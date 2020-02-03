
# (some of the test cases demonstrate the equations well)

# RPM normalization will remove any spike-in reads before generating RPM NFs

# for batch normalization, sample names must contain "rep" followed by the
#   replicate indication; if they're not named like this, users can use the
#   sample_names argument to make them conform



#' Calculating spike-in normalization factors
#'
#' Use \code{getSpikeInNFs} to obtain the spike-in normalization factors, or
#' \code{spikeInNormGRanges} to return the input GRanges objects with their
#' readcounts spike-in normalized.
#'
#' @param dataset.gr A GRanges object, or (more typically) a list of GRanges
#'   objects.
#' @param si_pattern A regular expression that matches spike-in chromosomes. Can
#'   be used in addition to, or as an alternative to \code{si_names}.
#' @param si_names A character vector giving the names of the spike-in
#'   chromosomes. Can be used in addition to, or as an alternative to
#'   \code{si_pattern}.
#' @param method One of the shown methods, which generate normalization factors
#'   for converting raw readcounts into "Spike-in normalized Reads Per Million
#'   mapped in Control" (the default), "Spike-in Normalized Read counts", or
#'   "Reads Per Million mapped". See descriptions below.
#' @param batch_norm A logical indicating if batch normalization should be used
#'   (\code{TRUE} by default). See descriptions below. If batch normalization is
#'   used, sample names must end with "rep#", wherein "#" is one or more
#'   characters (usually a number) giving the replicate. If this is not the
#'   case, users can use the \code{sample_names} argument to make the names
#'   conform.
#' @param ctrl_pattern A regular expression that matches negative control sample
#'   names.
#' @param ctrl_names A character vector giving the names of the negative control
#'   samples. Can be used as an alternative to \code{ctrl_pattern}.
#' @param field The metadata field in \code{dataset.gr} that contains raw
#'   readcounts. If each range is an individual read, set \code{field = NULL}.
#' @param sample_names An optional character vector that can be used to rename
#'   the samples in \code{dataset.gr}. Intended use is if \code{dataset.gr} is
#'   an unnamed list, or if \code{batch_norm = TRUE} but the sample names don't
#'   conform to the required naming scheme.
#' @param ncores The number of cores to use for computations.
#'
#' @return A numeric vector of normalization factors for each sample in
#'   \code{dataset.gr}. Normalization factors are to be applied by
#'   multiplication.
#'
#' @section Spike-in normalized Reads Per Million mapped in Control (SRPMC):
#'   This is the default spike-in normalization method, as its meaning is the
#'   most portable and generalizable. Experimental Reads Per Spike-in read (RPS)
#'   are calculated for each sample, \eqn{i}:
#'
#'   \deqn{RPS_i=\frac{experimental\_reads_i}{
#'   spikein\_reads_i}}{RPS_i = experimental.reads_i / spikein.reads_i}
#'
#'   RPS for each sample is divided by RPS for the negative control, which
#'   measures the change in total material vs. the negative control. This global
#'   adjustment is applied to standard RPM normalization for each sample:
#'
#'   \deqn{NF_i=\frac{RPS_i}{RPS_{control}} \cdot \frac{1 x
#'   10^6}{experimental\_reads_i}}{
#'   NF_i = (RPS_i / RPS_ctrl) x (1e6 / experimental.reads_i)}
#'
#'   Thus, the negative control(s) are simply RPM-normalized, while the other
#'   conditions are in equivalent, directly-comparable units ("Reads Per Million
#'   mapped reads in a negative control").
#'
#'   If \code{batch_norm = TRUE} (the default), all negative controls will be
#'   RPM-normalized, and the global changes in material for all other samples
#'   are calculated \emph{within each batch} (vs. the negative control within
#'   the same batch).
#'
#'   If \code{batch_norm = FALSE}, all samples are compared to the average RPS
#'   of the negative controls. This method can only be justified if batch has
#'   less effect on RPS than other sources of variation.
#'
#' @section Spike-in Normalized Reads (SNR): If \code{batch_norm = FALSE}, these
#'   normalization factors act to scale down the readcounts in each sample to
#'   make the spike-in read counts match the sample with the lowest number of
#'   spike-in reads:
#'
#'   \deqn{NF_i=\frac{min(spikein\_reads)}{spikein\_reads_i}}{
#'   NF_i = min(spikein.reads) / spikein.reads_i}
#'
#'   If \code{batch_norm = TRUE}, such normalization factors are calculated
#'   within each batch, but a final batch (replicate) adjustment is performed
#'   that results in the negative controls having the same normalized
#'   readcounts. In this way, the negative controls are used to adjust the
#'   normalized readcounts of their entire replicate. Just as when
#'   \code{batch_norm = FALSE}, one of the normalization factors will be
#'   \code{1}, while the rest will be \code{<1}.
#'
#'   One use for these normalization factors is for normalizing-by-subsampling;
#'   see \code{\link[BRGenomics:subsampleBySpikeIn]{subsampleBySpikeIn}}.
#'
#' @section Reads Per Million mapped reads (RPM): A simple convenience wrapper
#'   for calculating normalization factors for RPM normalization:
#'
#'   \deqn{NF_i=\frac{1 x 10^6}{experimental\_reads_i}}{
#'   NF_i = 1e6 / experimental.reads_i}
#'
#'   If spike-in reads are present, they're removed before the normalization
#'   factors are calculated.
#'
#' @author Mike DeBerardine
#' @seealso \code{\link[BRGenomics:getSpikeInCounts]{getSpikeInCounts}}
#' @export
#' @importFrom parallel detectCores
#'
#' @examples
getSpikeInNFs <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                          method = c("SRPMC", "SNR", "RPM"), batch_norm = TRUE,
                          ctrl_pattern = NULL, ctrl_names = NULL,
                          field = "score", sample_names = NULL,
                          ncores = detectCores()) {

    method <- match.arg(method, c("SRPMC", "SNR", "RPM"))

    if (!is.null(sample_names)) {
        if (!is.list(dataset.gr))  dataset.gr <- list(dataset.gr)
        names(dataset.gr) <- sample_names
    }

    counts.df <- getSpikeInCounts(dataset.gr, si_pattern, si_names, field,
                                  ncores)

    if (method == "SRPMC") {
        .get_nf_srpmc(counts.df, ctrl_pattern, ctrl_names, batch_norm)
    } else if (method == "SNR") {
        .get_nf_snr(counts.df, ctrl_pattern, ctrl_names, batch_norm)
    } else {
        1e6 / counts.df$exp_reads # RPM normalize
    }
}



.get_nf_snr <- function(df, ctrl_pattern, ctrl_names, batch_norm) {

    if (!batch_norm) {
        return( min(df$spike_reads) / df$spike_reads )

    } else {
        srep <- sub(".*rep", "rep", df$sample)

        # 1. get spike-in NFs within each replicate;
        # 2. normalize the replicates to one another, based on the controls
        snf <- rep(NA, nrow(df)) # initialize
        repnf <- rep(NA, nrow(df))

        idx_ctrl <- .get_idx_ctrl(df, ctrl_pattern, ctrl_names)
        min_ctrl <- idx_ctrl[ which.min(df$exp_reads[idx_ctrl]) ]
        min_rep <- srep[min_ctrl] # batch with lowest control reads

        for (i in unique(srep)) {
            idx_i <- which(srep == i)
            idx_ci <- intersect(idx_i, idx_ctrl)
            spike_i <- df$spike_reads[idx_i]

            # get within-replicate spike-in NF
            snf[idx_i] <- min(spike_i) / spike_i

            # get across-replicate NF (based on negative controls)
            repnf[idx_i] <- df$exp_reads[min_ctrl] / df$exp_reads[idx_ci]
        }

        # multiply within-batch NF and across-batch NFs
        combnf <- snf * repnf

        # adjust NFs in case all are now less than 1
        return(combnf / max(combnf))
    }
}


.get_nf_srpmc <- function(counts.df, ctrl_pattern, ctrl_names, batch_norm) {
    idx_ctrl <- .get_idx_ctrl(counts.df, ctrl_pattern, ctrl_names)
    rps <- counts.df$exp_reads / counts.df$spike_reads

    if (batch_norm) {
        srep <- sub(".*rep", "rep", counts.df$sample)
        ratio_ctrl <- rep(NA, nrow(counts.df))
        for (i in unique(srep)) {
            idx_i <- which(srep == i)
            idx_ci <- intersect(idx_i, idx_ctrl)
            ratio_ctrl[idx_i] <- rps[idx_i] / rps[idx_ci]
        }
    } else {
        ctrl_avg <- mean(rps[idx_ctrl])
        ratio_ctrl <- rps / ctrl_avg
    }

    nf_rpm <- 1e6 / counts.df$exp_reads
    return(ratio_ctrl * nf_rpm)
}


.get_idx_ctrl <- function(counts.df, ctrl_pattern, ctrl_names) {
    .check_xor_args(ctrl_pattern, ctrl_names)
    if (is.null(ctrl_names))
        ctrl_names <- grep(ctrl_pattern, counts.df$sample, value = TRUE)
    which(counts.df$sample %in% ctrl_names)
}


#' @importFrom parallel detectCores mcMap
#' @rdname getSpikeInNFs
spikeInNormGRanges <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                               method = c("SRPMC", "SNR", "RPM"),
                               batch_norm = TRUE, ctrl_pattern = NULL,
                               ctrl_names = NULL, field = "score",
                               sample_names = NULL, ncores = detectCores()) {

    nf <- getSpikeInNFs(dataset.gr, si_pattern, si_names, method, batch_norm,
                        ctrl_pattern, ctrl_names, field, sample_names, ncores)

    if (is.null(field)) {
        warning("With field = NULL, will calculate stranded coverage and make
                GRanges basepair-resolution before returning normalized
                GRanges")
        dataset.gr <- mclapply(dataset.gr, function(x) {
            makeGRangesBRG(getStrandedCoverage(x, field = NULL))
        }, mc.cores = ncores)
    }

    mcMap(function(x, field, nf) {
        gr <- x
        mcols(gr)[field] <- mcols(gr)[[field]] * nf
        gr
    }, dataset.gr, field, nf, mc.cores = ncores)
}

# goal of SNR (spike-in normalized reads) is to downward adjust readcounts to be
# correctly normalized; and then this can be used to subsample reads
#   AFTER doing the subsampling, however, we could adjust the units to be in RPM;
#   since all samples are then 1:1 in terms of raw (subsampled) readcounts,
#   we would just use the negative control's read-depth to RPM normalize all
#   samples

# dataset.gr should be a list; but but should add support for using multiplex
#   (already have a helper function somewhere for that...)

# RPM_units will convert the raw readcounts into RPM, using the negative control
#   if batch normalization is used, all negative controls have the same number
#   of adjusted readcounts, and that number is used;
#   if batch_norm = FALSE, the average number of reads across all negative
#   controls is used


#' Randomly subsample reads according to spike-in normalization
#'
#' @param
#' dataset.gr,si_pattern,si_names,ctrl_pattern,ctrl_names,batch_norm,field,sample_names,ncores
#' See \code{\link[BRGenomics:getSpikeInNFs]{getSpikeInNFs}}
#' @param RPM_units If set to \code{TRUE}, the final readcount values will be
#'   converted to \code{negative control RPM}.
#'
#' @return A list of GRanges parallel to \code{dataset.gr}, but containing fewer
#'   reads.
#' @author Mike DeBerardine
#'
#' @seealso \code{\link[BRGenomics:getSpikeInCounts]{getSpikeInCounts}},
#'   \code{\link[BRGenomics:getSpikeInNFs]{getSpikeInNFs}}
#'
#' @importFrom parallel detectCores mcMap
#' @export
#'
#' @examples
subsampleBySpikeIn <- function(dataset.gr, si_pattern = NULL, si_names = NULL,
                               ctrl_pattern = NULL, ctrl_names = NULL,
                               batch_norm = TRUE, RPM_units = FALSE,
                               field = "score", sample_names = NULL,
                               ncores = detectCores()) {

    if (!is.null(sample_names)) {
        if (!is.list(dataset.gr))  dataset.gr <- list(dataset.gr)
        names(dataset.gr) <- sample_names
    }

    counts.df <- getSpikeInCounts(dataset.gr, si_pattern, si_names, field,
                                  ncores)

    # get NFs for subsampling
    nf_snr <- .get_nf_snr(counts.df, ctrl_pattern, ctrl_names, batch_norm)
    nreads <- floor(counts.df$exp_reads * nf_snr)

    samples.gr <- mcMap(subsampleGRanges, dataset.gr, n = nreads, field = field,
                        mc.cores = ncores)

    if (RPM_units) {
        # get RPM NF for negative controls
        idx_ctrl <- .get_idx_ctrl(counts.df, ctrl_pattern, ctrl_names)
        nf_rpm <- 1e6 / mean(nreads[idx_ctrl])
        samples.gr <- mcMap(function(x, field) {
            gr <- x
            mcols(gr)[field] <- mcols(gr)[[field]] * nf_rpm
            gr
        }, samples.gr, field, mc.cores = ncores)
    }

    return(samples.gr)
}
