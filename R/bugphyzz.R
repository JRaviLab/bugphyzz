utils::globalVariables(c(
    "Rank",
    "Attribute", "Attribute_group", "Attribute_new", "Attribute_range",
    "Attribute_value", "Attribute_value_max", "Attribute_value_min", "Br-C10:1",
    "Evidence", "Frequency", "NCBI_ID", "Oxo-C19:1", "Taxon_name", "Unit",
    "attribute", "functionname", "gram_stain", "physiology", "unit", "value"
))

#' Import bugphyzz
#'
#' \code{importBugphyzz} imports bugphyzz annotations as a list of
#' tidy data.frames. To learn more about the structure of the data.frames
#' please check the bugphyzz vignette with `browseVignettes("bugphyzz")`.
#'
#' @param version Character string indicating the version. Default is the
#' latest release on Zenodo.  Options: Zenodo DOI, GitHub commit hash, or devel.
#' @param forceDownload Logical value. Force a fresh download of the data or
#' use the one stored in the cache (if available). Default is FALSE.
#' @param v Validation value. Default 0.8 (see details).
#' @param excludeRarely Default is TRUE. Exclude values with
#' Frequency == FALSE (see details).
#'
#' @details
#'
#' ## Data structure
#' The data structure of the data.frames imported with `importBugphyzz` are
#' detailed in the main vignette. Please run `browseVignettes("bugphyzz")`.
#'
#' ## Validation (`v` argument)
#' Data imported with `importBugphyzz` includes annotations imputed through
#' ancestral state reconstruction (ASR) methods. A 10-fold cross-validation
#' approach was implemented to assess the reliability of the data imputed.
#' Mathew's correlation coefficient (MCC) and R-squared (R2) were used for the
#' validation of discrete and numeric attributes.
#' Details can be found at: https://github.com/waldronlab/taxPProValidation.
#' By default, imputed annotations with a MCC or R2 value greater than 0.5 are
#' imported. The minimum value can be adjusted with the `v` argument (only
#' values between 0 and 1).
#'
#' ## Frequency (excludeRarely argument)
#' One of the variables in the bugphyzz data.frames is "Frequency", which
#' can adopt values of
#' "always", "usually", "sometimes", "rarely", or "never". By default
#' "never" and "rarely" are excluded. "rarely" could be included with
#' `excludeRarely = FALSE`. To learn more about these frequency keywords
#' please check the bugphyzz vignette with `browseVignettes("bugphyzz")`.
#'
#' @return A list of tidy data frames.
#' @export
#'
#' @examples
#'
#' bp <- importBugphyzz()
#' names(bp)
#'
importBugphyzz <- function(
        version = "10.5281/zenodo.10980813", forceDownload = FALSE, v = 0.8,
        excludeRarely = TRUE

) {

    ## output is a list of three data.frames
    ## one of each: binary, multistate, numeric
    output <- .downloadResource(version, forceDownload)

    ## TODO add release version
    output <- lapply(output, function(x) split(x, x$Attribute))
    output <- purrr::list_flatten(output)

    ## TODO correct plant pathogenicity name earlier in the workflow or
    ## better yet, directly in the curation
    pos <- which(names(output) == "plant pathogenity")
    names(output)[pos] <- "plant pathogenicity"
    output <- purrr::map(output, ~ {
        .x |>
            dplyr::mutate(
                Attribute = ifelse(
                    Attribute == "plant pathogenity",
                    "plant pathogenicity",
                    Attribute
                )
            )
    })

    names(output) <- purrr::map_chr(output, ~ unique(.x$Attribute))
    val <- .validationData() |>
        dplyr::filter(rank == "all") |>
        dplyr::select(physiology, attribute, value) |>
        dplyr::mutate(physiology = tolower(physiology)) |>
        dplyr::mutate(attribute = tolower(attribute))

    output <- purrr::map(output, ~ {
        attrType <- unique(.x$Attribute_type)
        if (attrType == "binary") {
            val <- dplyr::select(val, Attribute = attribute, value)
            o <- dplyr::left_join(.x, val, by = "Attribute" )
        } else if (
            attrType == "multistate-intersection" ||
            attrType == "multistate-union"
        ) {
            val <- dplyr::select(
                val, Attribute = physiology, Attribute_value = attribute, value
            )
            o <- dplyr::left_join(
                dplyr::mutate(.x, Attribute_value = tolower(Attribute_value)),
                val, by = c("Attribute", "Attribute_value")
            )
        } else if (attrType == "numeric") {
            val <- dplyr::select(val, Attribute = attribute, value)
            o <- dplyr::left_join(.x, val, by = "Attribute") |>
                dplyr::rename(NSTI = nsti)
        }
        o |>
            dplyr::filter(
                !(value < v & Evidence == "asr")
            ) |>
            dplyr::mutate(value = ifelse(Evidence != "asr", NA, value)) |>
            dplyr::rename(Validation = value)
    })

    if (excludeRarely) {
        output <- purrr::map(
            output, ~ dplyr::filter(.x, Frequency != "rarely")
        )
    }
    return(output)
}

#' Make signatures
#'
#' \code{makeSignatures} Creates signatures for a list of bug signatures from
#' a tidy data.frame imported through the `importBugphyzz` function. Please
#' run `browseVignettes("bugphyz")` for detailed examples.
#'
#' @param dat A data.frame.
#' @param taxIdType A character string. Valid options: NCBI_ID, Taxon_name.
#' @param taxLevel A character vector. Taxonomic rank. Valid options:
#' superkingdom, kingdom, phylum, class, order, family, genus, species, strain.
#' They can be combined. "mixed" is equivalent to select all valid ranks.
#' @param evidence A character vector. Valid options: exp, igc, nas, tas, tax,
#' asr. They can be combined. Default is all.
#' @param frequency A character vector. Valid options: always, usually,
#' sometimes, rarely, unknown. They can be combined. By default, "rarely" is
#' excluded.
#' @param minSize Minimum number of bugs in a signature. Default is 10.
#' @param min Minimum value (inclusive). Only for numeric attributes.
#' Default is NULL.
#' @param max Maximum value (inclusive). Only for numeric attributes.
#' Default is NULL.
#'
#' @return A list of character vectors with scientific names or taxids.
#' @export
#'
#' @examples
#'
#' bp <- importBugphyzz()
#' sigs <- purrr::map(bp, makeSignatures)
#' sigs <- purrr::list_flatten(sigs, name_spec = "{inner}")
#'
makeSignatures <- function(
        dat, taxIdType = "NCBI_ID",
        taxLevel = "mixed",
        evidence = c("exp", "igc", "tas", "nas", "tax", "asr"),
        frequency = c("always", "usually", "sometimes", "unknown"),
        minSize = 10, min = NULL, max = NULL
) {
    attrType <- unique(dat$Attribute_type)
    if ("mixed" %in% taxLevel) {
        taxLevel <- c(
            "kingdom", "phylum", "class", "order", "family", "genus", "species",
            "strain"
        )
    }
    dat <- dat |>
        {\(y) y[which(y$Rank %in% taxLevel),]}() |> 
        {\(y) y[which(y$Evidence %in% evidence),]}() |> 
        {\(y) y[which(y$Frequency %in% frequency),]}()
        
    if (!nrow(dat)) {
        warning(
            "Not enough data for creating signatures.",
            " Try different filtering options",
            call. = FALSE
        )
        return(NULL)
    }
    if (
        attrType %in%
        c("multistate-intersection", "binary", "multistate-union")
        ) {
        s <- .makeSignaturesDiscrete(dat = dat, taxIdType = taxIdType)
    } else if (attrType %in% c("range", "numeric")) {
        s <- .makeSignaturesNumeric(
            dat = dat, taxIdType = taxIdType, min = min, max = max
        )
    }
    output <- purrr::keep(s, ~ length(.x) >= minSize)
    if (!length(output)) {
        warning(
            "Not enough data for creating signatures.",
            " Try different filtering options",
            call. = FALSE
        )
    }
    return(output)
}

#' Get Taxon Signatures
#'
#' \code{getTaxonSignatures} returns the names of all of the signatures
#' associated with a particular taxon. More details can be found in the main
#' bugphyzz vignette; please run `browseVignettes("bugphyzz")`.
#'
#' @param tax A valid NCBI ID or taxon name. If taxon name is used, the
#' argument taxIdType = "Taxon_name" must also be used.
#' @param bp List of data.frames imported with \code{importBugphyzz}.
#' @param ... Arguments passed to \code{makeSignatures}.
#'
#' @return A character vector with the names of the signatures for a taxon.
#' @export
#'
#' @examples
#' taxid <- "562"
#' taxonName <- "Escherichia coli"
#' bp <- importBugphyzz()
#' sig_names_1 <- getTaxonSignatures(taxid, bp)
#' sig_names_2 <- getTaxonSignatures(taxonName, bp, taxIdType = "Taxon_name")
#'
getTaxonSignatures <- function(tax, bp, ...) {
    sigs <- purrr::map(bp, makeSignatures, ...)
    sigs <- purrr::list_flatten(sigs, name_spec = "{inner}")
    pos <- which(purrr::map_lgl(sigs, ~ tax %in% .x))
    output <- names(sigs)[pos]
    return(output)
}

# Non exported functions ----------------------------------------------------
.makeSignaturesDiscrete <- function(dat, taxIdType = "NCBI_ID") {
    dat$Attribute <- paste0(
        "bugphyz:", dat$Attribute, "|", dat$Attribute_value
    )
    dat |> 
        {\(y) S4Vectors::split(y, y$Attribute)}() |>
        lapply(function(x) unique(x[[taxIdType]]))
}

.makeSignaturesNumeric <- function(
        dat, taxIdType = "NCBI_ID", min = NULL, max = NULL
) {
    if (!is.null(min) || !is.null(max)) {
        if (is.null(min)) {
            message(
                "Minimum unespecified. Using ", min(dat$Attribute_value), "."
            )
            min <- min(dat$Attribute_value)
        }
        if (is.null(max)) {
            message(
                "Maximum unespecified. Using ", max(dat$Attribute_value), "."
            )
            max <- max(dat$Attribute_value)
        }
        
        dat <- dat[
            which(dat$Attribute_value >= min & dat$Attribute_value <= max),
        ]
        dat$Attribute <- paste0(
            "bugphyzz:", dat$Attribute, "| >=", min, " & <=", max
        )
    } else {
        thr <- .thresholds() |>
            dplyr::filter(Attribute_group == unique(dat$Attribute))
        attrName <- thr$Attribute
        minValues <- thr$lower
        maxValues <- thr$upper
        dat$tmp_col <- NA
        for (i in seq_along(attrName)) {
            if (is.na(minValues[i]))
                minValues[i] <- min(dat$Attribute_value) - 0.01
            if (is.na(maxValues[i]))
                maxValues[i] <- max(dat$Attribute_value)
            pos <- which(
                dat$Attribute_value > minValues[i] &
                    dat$Attribute_value <= maxValues[i]
            )
            dat$tmp_col[pos] <- attrName[i]
            dat$Attribute[pos] <- paste0(
                "bugphyzz:", dat$Attribute[pos], "|", attrName[i], "| > ",
                round(minValues[i], 2), " & <= ", maxValues[i]
            )
        }
    }
    dat |>
        {\(y) S4Vectors::split(y, y$Attribute)}() |>
        lapply(function(x) unique(x[[taxIdType]]))
}

.thresholds <- function() {
    fpath <- file.path('extdata', 'thresholds.tsv')
    fname <- system.file(fpath, package = 'bugphyzz', mustWork = TRUE)
    utils::read.table(fname, header = TRUE, sep = '\t') |>
        dplyr::mutate(
            range = dplyr::case_when(
                is.na(lower) ~ paste0('<=', upper),
                is.na(upper) ~ paste0('>=', lower),
                TRUE ~ paste0(lower, '-', upper)
            ),
            unit = ifelse(is.na(unit), '', unit)
        ) |>
        dplyr::mutate(Attribute_range = paste0(range, unit)) |>
        dplyr::relocate(
            Attribute_group, Attribute, Attribute_range
        )
}

.validationData <- function() {
    fname <- system.file(
        "extdata", "validation_summary.tsv", package = "bugphyzz"
    )
    utils::read.table(
        file = fname, header = TRUE, sep = "\t", row.names = NULL
    ) |>
        dplyr::mutate(
            value = dplyr::case_when(
                !is.na(mcc_mean) & is.na(r2_mean) ~ mcc_mean,
                is.na(mcc_mean) & !is.na(r2_mean) ~ r2_mean
            )
        )
}

## Import a version of bupghyzz
.downloadResource <- function(version, forceDownload) {
    if (stringr::str_detect(version, "^10.5281/zenodo.[0-9]+$")) {
        suffix <- sub("^10.5281/zenodo\\.", "", version)
        output <- .downloadZ(suffix, forceDownload)
    } else if (
        version == "devel" ||
            stringr::str_detect(version, stringr::regex("^[:alnum:]{7}$"))
        ){
        output <- .downloadGH(version, forceDownload)
    } else {
        stop("Version must be a Zenodo DOI, GitHub commit hash, or 'devel'.")
    }
    return(output)
}

## Function for downloading data on Zenodo
.downloadZ <- function(record, forceDownload) {
    baseUrl <- paste0("https://zenodo.org/api/records/", record)
    req <- httr2::request(baseUrl)
    res <- httr2::req_perform(req)
    l <- httr2::resp_body_json(res)

    fileNamesApi <- purrr::map_chr(l$files, ~ .x$links$self)
    fileNamesUrl <- sub(
        "(^.*)(api/)(.*)(/content$)", "\\1\\3", fileNamesApi
    )

    rpath <- .getResource(
        rname = paste0("bugphyzz.zip"),
        url = fileNamesUrl, verbose = TRUE, force = forceDownload
    )
    tempDir <- tempdir()
    utils::unzip(zipfile = rpath, exdir = tempDir, junkpaths = TRUE)
    files <- list.files(tempDir, pattern = "csv", full.names = TRUE)

    output <- vector("list", length(files))
    for (i in seq_along(output)) {
        output[[i]] <- utils::read.csv(files[i], header = TRUE, skip = 1) |>
            dplyr::mutate(Attribute = tolower(Attribute))
    }
    return(output)
}

## Function for downloading data on GitHub
.downloadGH <- function(version, forceDownload) {
    fileSuffix <- c("binary", "multistate", "numeric")
    urls <- paste0(
        "https://github.com/waldronlab/bugphyzzExports/raw/", version,
        "/bugphyzz_", fileSuffix, ".csv"
    )
    names(urls) <-  c("binary", "multistate", "numeric")
    output <- vector("list", length(urls))
    for (i in seq_along(output)) {
        message("Importing ", names(urls)[i], " data...")
        names(output)[i] <- names(urls)[i]
        rpath <- .getResource(
            rname = paste0("bugphyzz_", names(urls)[i], ".csv"),
            url = urls[i], verbose = TRUE, force = forceDownload
        )
        output[[i]] <- utils::read.csv(rpath, header = TRUE, skip = 1) |>
            dplyr::mutate(Attribute = tolower(Attribute))
    }
    return(output)
}
