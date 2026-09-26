# ------------------------------------------------------------
# Purpose:
# Perform pooled functional enrichment analysis comparing
# MIsoTEx Ubiquitous Uniform vs Non-uniform transcripts.
#
# All five expression types (Very High, High, Moderate, Low,
# Very Low) are pooled within Uniform and Non-uniform.
#
# Analyses:
# GO-BP, GO-MF, GO-CC, KEGG and Reactome
#
# Outputs:
# Enrichment tables, significant tables, transcript mappings,
# Uniform vs Non-uniform comparison tables, ellipse-shaped
# Venn diagrams, summary tables and final summary image.
#
# Venn:
# Uniform / UUni     = Pink
# Non-uniform / UVar = Grey
#
# No >50 transcript filtering is applied.
#
# Author: SBM
# ------------------------------------------------------------


# ============================================================
# 1. USER SETTINGS
# ============================================================

UNIFORM_FILE <-
    "Uniform_quantity_of_ubiquitously_present_list_Expression_profile_Score.csv"

NONUNIFORM_FILE <-
    "Varying_quantity_of_ubiquitously_present_list_Expression_profile_Score.csv"

OUTPUT_DIR <-
    "MIsoTEx_Ubiquitous_Uniform_vs_NonUniform_Pooled_Functional_Analysis"

FDR_CUTOFF <- 0.05


# ============================================================
# 2. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
    "clusterProfiler",
    "enrichplot",
    "KEGGREST",
    "ReactomePA",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "venn",
    "ggplot2",
    "dplyr",
    "grid"
)


cat("\n")
cat("============================================================\n")
cat("CHECKING REQUIRED PACKAGES\n")
cat("============================================================\n\n")


missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]


if (length(missing_packages) > 0) {

    cat("MISSING PACKAGES:\n\n")
    print(missing_packages)

    stop(
        "\nPlease install the missing packages first.\n",
        "For ellipse Venn diagrams use:\n",
        "install.packages(\"venn\")\n"
    )
}


cat("All required packages are installed.\n")


# ============================================================
# 3. LOAD PACKAGES
# ============================================================

suppressPackageStartupMessages({

    library(clusterProfiler)
    library(enrichplot)
    library(KEGGREST)
    library(ReactomePA)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(venn)
    library(ggplot2)
    library(dplyr)
    library(grid)

})


# ============================================================
# 4. PACKAGE VERSIONS
# ============================================================

cat("\n")
cat("PACKAGE VERSIONS\n")
cat("----------------\n")


for (pkg in required_packages) {

    cat(
        sprintf(
            "%-20s : %s\n",
            pkg,
            as.character(
                packageVersion(pkg)
            )
        )
    )

}


# ============================================================
# 5. CHECK INPUT FILES
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CHECKING INPUT FILES\n")
cat("============================================================\n")


if (!file.exists(UNIFORM_FILE)) {

    stop(
        "\nUniform file not found:\n",
        UNIFORM_FILE,
        "\n"
    )

}


if (!file.exists(NONUNIFORM_FILE)) {

    stop(
        "\nNon-uniform file not found:\n",
        NONUNIFORM_FILE,
        "\n"
    )

}


cat("\nUniform file:\n")
cat(UNIFORM_FILE, "\n")


cat("\nNon-uniform file:\n")
cat(NONUNIFORM_FILE, "\n")


# ============================================================
# 6. CREATE OUTPUT DIRECTORIES
# ============================================================

dirs <- c(

    "",

    "Uniform",

    "Uniform/GO_BP",
    "Uniform/GO_MF",
    "Uniform/GO_CC",
    "Uniform/KEGG",
    "Uniform/Reactome",

    "Non_Uniform",

    "Non_Uniform/GO_BP",
    "Non_Uniform/GO_MF",
    "Non_Uniform/GO_CC",
    "Non_Uniform/KEGG",
    "Non_Uniform/Reactome",

    "Uniform_vs_NonUniform",

    "Uniform_vs_NonUniform/GO_BP",
    "Uniform_vs_NonUniform/GO_MF",
    "Uniform_vs_NonUniform/GO_CC",
    "Uniform_vs_NonUniform/KEGG",
    "Uniform_vs_NonUniform/Reactome",

    "Plots"

)


for (d in dirs) {

    dir.create(
        file.path(
            OUTPUT_DIR,
            d
        ),
        recursive = TRUE,
        showWarnings = FALSE
    )

}


# ============================================================
# 7. READ INPUT
# ============================================================

read_input <- function(
    file,
    group_name
) {

    cat("\n")
    cat("------------------------------------------------------------\n")
    cat(
        "READING ",
        group_name,
        "\n",
        sep = ""
    )
    cat("------------------------------------------------------------\n")


    df <- read.csv(
        file,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )


    cat(
        "Total rows: ",
        nrow(df),
        "\n"
    )


    cat("Columns:\n")
    print(colnames(df))


    # --------------------------------------------------------
    # Standardize column names
    # --------------------------------------------------------

    if ("Gene_stable_ID" %in% colnames(df)) {

        colnames(df)[
            colnames(df) == "Gene_stable_ID"
        ] <- "gene_stable_id"

    }


    if ("Gene_name" %in% colnames(df)) {

        colnames(df)[
            colnames(df) == "Gene_name"
        ] <- "gene_name"

    }


    # --------------------------------------------------------
    # Required columns
    # --------------------------------------------------------

    required_columns <- c(
        "transcript_id",
        "gene_stable_id",
        "gene_name",
        "Expression_type"
    )


    missing_columns <- setdiff(
        required_columns,
        colnames(df)
    )


    if (length(missing_columns) > 0) {

        stop(
            "\nMissing required columns in ",
            group_name,
            ":\n",
            paste(
                missing_columns,
                collapse = ", "
            ),
            "\n"
        )

    }


    # --------------------------------------------------------
    # Remove transcript version
    # --------------------------------------------------------

    df$transcript_id <- sub(
        "\\..*$",
        "",
        as.character(
            df$transcript_id
        )
    )


    # --------------------------------------------------------
    # Remove Ensembl gene version
    # --------------------------------------------------------

    df$gene_stable_id <- sub(
        "\\..*$",
        "",
        as.character(
            df$gene_stable_id
        )
    )


    # --------------------------------------------------------
    # Clean expression type
    # --------------------------------------------------------

    df$Expression_type <- trimws(
        as.character(
            df$Expression_type
        )
    )


    df$Group <- group_name


    return(df)

}


# ============================================================
# 8. READ BOTH FILES
# ============================================================

uniform_df <- read_input(
    UNIFORM_FILE,
    "Uniform"
)


nonuniform_df <- read_input(
    NONUNIFORM_FILE,
    "Non_Uniform"
)


# ============================================================
# 9. REMOVE DUPLICATE TRANSCRIPTS
# ============================================================

uniform_df <- uniform_df %>%

    distinct(
        transcript_id,
        .keep_all = TRUE
    )


nonuniform_df <- nonuniform_df %>%

    distinct(
        transcript_id,
        .keep_all = TRUE
    )


cat("\n")
cat("Uniform transcripts: ")
cat(nrow(uniform_df))
cat("\n")


cat("Non-uniform transcripts: ")
cat(nrow(nonuniform_df))
cat("\n")


# ============================================================
# 10. EXPRESSION TYPE DISTRIBUTION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("EXPRESSION TYPE DISTRIBUTION\n")
cat("============================================================\n")


cat("\nUNIFORM:\n")

print(
    table(
        uniform_df$Expression_type
    )
)


cat("\nNON-UNIFORM:\n")

print(
    table(
        nonuniform_df$Expression_type
    )
)


# ============================================================
# 11. SAVE CLEAN INPUTS
# ============================================================

write.csv(
    uniform_df,
    file.path(
        OUTPUT_DIR,
        "Uniform",
        "Uniform_All_Expression_Types.csv"
    ),
    row.names = FALSE
)


write.csv(
    nonuniform_df,
    file.path(
        OUTPUT_DIR,
        "Non_Uniform",
        "NonUniform_All_Expression_Types.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 12. COMBINE ONLY FOR MAPPING
# ============================================================

all_gene_ids <- unique(
    c(
        uniform_df$gene_stable_id,
        nonuniform_df$gene_stable_id
    )
)


cat("\n")
cat(
    "Total unique Ensembl genes across both groups: "
)

cat(
    length(all_gene_ids)
)

cat("\n")


# ============================================================
# 13. ENSEMBL → ENTREZ
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ENSEMBL → ENTREZ MAPPING\n")
cat("============================================================\n")


gene_map <- clusterProfiler::bitr(
    all_gene_ids,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)


gene_map <- gene_map %>%

    distinct(
        ENSEMBL,
        .keep_all = TRUE
    )


cat(
    "Mapped Ensembl genes: ",
    nrow(gene_map),
    "\n"
)


cat(
    "Mapping percentage: ",
    round(
        100 *
            nrow(gene_map) /
            length(all_gene_ids),
        2
    ),
    "%\n"
)


write.csv(
    gene_map,
    file.path(
        OUTPUT_DIR,
        "Uniform_vs_NonUniform",
        "Ensembl_to_Entrez_Mapping.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 14. ADD ENTREZ IDS
# ============================================================

uniform_map <- uniform_df %>%

    left_join(
        gene_map,
        by = c(
            "gene_stable_id" = "ENSEMBL"
        )
    )


nonuniform_map <- nonuniform_df %>%

    left_join(
        gene_map,
        by = c(
            "gene_stable_id" = "ENSEMBL"
        )
    )


# ============================================================
# 15. BACKGROUND
# ============================================================

background_genes <- AnnotationDbi::keys(
    org.Hs.eg.db,
    keytype = "ENTREZID"
)


background_genes <- unique(
    background_genes
)


cat(
    "\nBackground genes: ",
    length(background_genes),
    "\n"
)


# ============================================================
# 16. GET POOLED GENES
# ============================================================

uniform_genes <- unique(
    uniform_map$ENTREZID[
        !is.na(
            uniform_map$ENTREZID
        )
    ]
)


nonuniform_genes <- unique(
    nonuniform_map$ENTREZID[
        !is.na(
            nonuniform_map$ENTREZID
        )
    ]
)


cat("\n")
cat("============================================================\n")
cat("POOLED GENE SETS\n")
cat("============================================================\n")


cat("\nUniform:\n")

cat(
    "  Transcripts = ",
    nrow(uniform_map),
    "\n"
)

cat(
    "  Unique genes = ",
    length(uniform_genes),
    "\n"
)


cat("\nNon-uniform:\n")

cat(
    "  Transcripts = ",
    nrow(nonuniform_map),
    "\n"
)

cat(
    "  Unique genes = ",
    length(nonuniform_genes),
    "\n"
)


# ============================================================
# 17. TRANSCRIPT MAPPING FUNCTION
# ============================================================

get_transcripts <- function(
    gene_string,
    df
) {

    if (
        is.na(gene_string) ||
        gene_string == ""
    ) {

        return("")

    }


    ids <- unlist(
        strsplit(
            gene_string,
            "/"
        )
    )


    tmp <- df %>%

        filter(
            ENTREZID %in% ids
        )


    paste(
        unique(
            tmp$transcript_id
        ),
        collapse = ";"
    )

}


# ============================================================
# 18. RUN GO ENRICHMENT
# ============================================================

run_GO <- function(
    genes,
    df,
    group_name,
    ontology
) {

    cat(
        "\nRunning GO-",
        ontology,
        " for ",
        group_name,
        "\n",
        sep = ""
    )


    result <- tryCatch(

        {

            enrichGO(

                gene =
                    genes,

                universe =
                    background_genes,

                OrgDb =
                    org.Hs.eg.db,

                keyType =
                    "ENTREZID",

                ont =
                    ontology,

                pAdjustMethod =
                    "BH",

                pvalueCutoff =
                    1,

                qvalueCutoff =
                    1,

                readable =
                    FALSE

            )

        },

        error = function(e) {

            warning(
                "GO error: ",
                e$message
            )

            NULL

        }

    )


    if (
        is.null(result) ||
        nrow(
            as.data.frame(result)
        ) == 0
    ) {

        return(NULL)

    }


    result <- as.data.frame(
        result
    )


    result$Group <-
        group_name


    result$FDR_BH <-
        result$p.adjust


    result$Significant <-
        ifelse(
            !is.na(
                result$FDR_BH
            ) &
                result$FDR_BH < FDR_CUTOFF,
            "Yes",
            "No"
        )


    result$Transcript_IDs <- sapply(
        result$geneID,
        get_transcripts,
        df = df
    )


    result$Transcript_Count <- sapply(
        result$Transcript_IDs,
        function(x) {

            if (
                is.na(x) ||
                x == ""
            ) {

                return(0)

            }


            length(
                unique(
                    unlist(
                        strsplit(
                            x,
                            ";"
                        )
                    )
                )
            )

        }
    )


    return(result)

}


# ============================================================
# 19. RUN ALL GO
# ============================================================

UNIFORM_GO_BP <- run_GO(
    uniform_genes,
    uniform_map,
    "Uniform",
    "BP"
)


UNIFORM_GO_MF <- run_GO(
    uniform_genes,
    uniform_map,
    "Uniform",
    "MF"
)


UNIFORM_GO_CC <- run_GO(
    uniform_genes,
    uniform_map,
    "Uniform",
    "CC"
)


NONUNIFORM_GO_BP <- run_GO(
    nonuniform_genes,
    nonuniform_map,
    "Non_Uniform",
    "BP"
)


NONUNIFORM_GO_MF <- run_GO(
    nonuniform_genes,
    nonuniform_map,
    "Non_Uniform",
    "MF"
)


NONUNIFORM_GO_CC <- run_GO(
    nonuniform_genes,
    nonuniform_map,
    "Non_Uniform",
    "CC"
)


# ============================================================
# 20. SAVE RESULTS
# ============================================================

save_result <- function(
    result,
    group,
    database
) {

    if (
        is.null(result) ||
        nrow(result) == 0
    ) {

        return()

    }


    group_dir <- ifelse(
        group == "Uniform",
        "Uniform",
        "Non_Uniform"
    )


    write.csv(
        result,
        file.path(
            OUTPUT_DIR,
            group_dir,
            paste0(
                database,
                "_ALL.csv"
            )
        ),
        row.names = FALSE
    )


    significant_result <- result %>%

        filter(
            !is.na(FDR_BH),
            FDR_BH < FDR_CUTOFF
        )


    write.csv(
        significant_result,
        file.path(
            OUTPUT_DIR,
            group_dir,
            paste0(
                database,
                "_SIGNIFICANT.csv"
            )
        ),
        row.names = FALSE
    )

}


save_result(
    UNIFORM_GO_BP,
    "Uniform",
    "GO_BP"
)


save_result(
    UNIFORM_GO_MF,
    "Uniform",
    "GO_MF"
)


save_result(
    UNIFORM_GO_CC,
    "Uniform",
    "GO_CC"
)


save_result(
    NONUNIFORM_GO_BP,
    "Non_Uniform",
    "GO_BP"
)


save_result(
    NONUNIFORM_GO_MF,
    "Non_Uniform",
    "GO_MF"
)


save_result(
    NONUNIFORM_GO_CC,
    "Non_Uniform",
    "GO_CC"
)


# ============================================================
# 21. KEGG
# ============================================================

cat("\n")
cat("============================================================\n")
cat("KEGG ENRICHMENT\n")
cat("============================================================\n")


kegg_link <- tryCatch(

    KEGGREST::keggLink(
        "pathway",
        "hsa"
    ),

    error = function(e) NULL

)


if (!is.null(kegg_link)) {


    kegg_term2gene <- data.frame(

        Term_ID =
            sub(
                "path:",
                "",
                as.character(
                    kegg_link
                )
            ),

        ENTREZID =
            sub(
                "hsa:",
                "",
                names(
                    kegg_link
                )
            ),

        stringsAsFactors =
            FALSE

    )


    kegg_names <- tryCatch(

        KEGGREST::keggList(
            "pathway",
            "hsa"
        ),

        error = function(e) NULL

    )


    if (!is.null(kegg_names)) {

        kegg_term2name <- data.frame(

            Term_ID =
                sub(
                    "path:",
                    "",
                    names(
                        kegg_names
                    )
                ),

            Description =
                sub(
                    " - Homo sapiens \\(human\\)$",
                    "",
                    as.character(
                        kegg_names
                    )
                ),

            stringsAsFactors =
                FALSE

        )

    } else {

        kegg_term2name <- NULL

    }


    run_KEGG <- function(
        genes,
        df,
        group
    ) {


        result <- tryCatch(

            {

                enricher(

                    gene =
                        genes,

                    universe =
                        background_genes,

                    TERM2GENE =
                        kegg_term2gene,

                    TERM2NAME =
                        kegg_term2name,

                    pAdjustMethod =
                        "BH",

                    pvalueCutoff =
                        1,

                    qvalueCutoff =
                        1

                )

            },

            error = function(e) {

                warning(
                    "KEGG error: ",
                    e$message
                )

                NULL

            }

        )


        if (
            is.null(result) ||
            nrow(
                as.data.frame(result)
            ) == 0
        ) {

            return(NULL)

        }


        result <- as.data.frame(
            result
        )


        result$Group <-
            group


        result$FDR_BH <-
            result$p.adjust


        result$Significant <-
            ifelse(
                !is.na(
                    result$FDR_BH
                ) &
                    result$FDR_BH < FDR_CUTOFF,
                "Yes",
                "No"
            )


        result$Transcript_IDs <- sapply(
            result$geneID,
            get_transcripts,
            df = df
        )


        result$Transcript_Count <- sapply(
            result$Transcript_IDs,
            function(x) {

                if (
                    is.na(x) ||
                    x == ""
                ) {

                    return(0)

                }


                length(
                    unique(
                        unlist(
                            strsplit(
                                x,
                                ";"
                            )
                        )
                    )
                )

            }
        )


        return(result)

    }


    UNIFORM_KEGG <- run_KEGG(
        uniform_genes,
        uniform_map,
        "Uniform"
    )


    NONUNIFORM_KEGG <- run_KEGG(
        nonuniform_genes,
        nonuniform_map,
        "Non_Uniform"
    )


    save_result(
        UNIFORM_KEGG,
        "Uniform",
        "KEGG"
    )


    save_result(
        NONUNIFORM_KEGG,
        "Non_Uniform",
        "KEGG"
    )


} else {

    UNIFORM_KEGG <- NULL

    NONUNIFORM_KEGG <- NULL

    cat(
        "KEGG retrieval failed.\n"
    )

}


# ============================================================
# 22. REACTOME
# ============================================================

cat("\n")
cat("============================================================\n")
cat("REACTOME ENRICHMENT\n")
cat("============================================================\n")


run_reactome <- function(
    genes,
    df,
    group
) {


    result <- tryCatch(

        {

            ReactomePA::enrichPathway(

                gene =
                    genes,

                universe =
                    background_genes,

                organism =
                    "human",

                pAdjustMethod =
                    "BH",

                pvalueCutoff =
                    1,

                qvalueCutoff =
                    1,

                readable =
                    FALSE

            )

        },

        error = function(e) {

            warning(
                "Reactome error: ",
                e$message
            )

            NULL

        }

    )


    if (
        is.null(result) ||
        nrow(
            as.data.frame(result)
        ) == 0
    ) {

        return(NULL)

    }


    result <- as.data.frame(
        result
    )


    result$Group <-
        group


    result$FDR_BH <-
        result$p.adjust


    result$Significant <-
        ifelse(
            !is.na(
                result$FDR_BH
            ) &
                result$FDR_BH < FDR_CUTOFF,
            "Yes",
            "No"
        )


    result$Transcript_IDs <- sapply(
        result$geneID,
        get_transcripts,
        df = df
    )


    result$Transcript_Count <- sapply(
        result$Transcript_IDs,
        function(x) {

            if (
                is.na(x) ||
                x == ""
            ) {

                return(0)

            }


            length(
                unique(
                    unlist(
                        strsplit(
                            x,
                            ";"
                        )
                    )
                )
            )

        }
    )


    return(result)

}


UNIFORM_REACTOME <- run_reactome(
    uniform_genes,
    uniform_map,
    "Uniform"
)


NONUNIFORM_REACTOME <- run_reactome(
    nonuniform_genes,
    nonuniform_map,
    "Non_Uniform"
)


save_result(
    UNIFORM_REACTOME,
    "Uniform",
    "Reactome"
)


save_result(
    NONUNIFORM_REACTOME,
    "Non_Uniform",
    "Reactome"
)


# ============================================================
# 23. UNIFORM vs NON-UNIFORM COMPARISON
# ============================================================

compare_functions <- function(
    uniform_result,
    nonuniform_result,
    database
) {


    cat("\n")
    cat("============================================================\n")
    cat(
        "COMPARING: ",
        database,
        "\n"
    )
    cat("============================================================\n")


    # --------------------------------------------------------
    # Empty protection
    # --------------------------------------------------------

    if (
        is.null(uniform_result)
    ) {

        uniform_result <- data.frame(
            ID = character(),
            Description = character(),
            pvalue = numeric(),
            p.adjust = numeric(),
            FDR_BH = numeric(),
            stringsAsFactors = FALSE
        )

    }


    if (
        is.null(nonuniform_result)
    ) {

        nonuniform_result <- data.frame(
            ID = character(),
            Description = character(),
            pvalue = numeric(),
            p.adjust = numeric(),
            FDR_BH = numeric(),
            stringsAsFactors = FALSE
        )

    }


    # --------------------------------------------------------
    # Significant terms
    # --------------------------------------------------------

    U <- uniform_result %>%

        filter(
            !is.na(FDR_BH),
            FDR_BH < FDR_CUTOFF
        )


    N <- nonuniform_result %>%

        filter(
            !is.na(FDR_BH),
            FDR_BH < FDR_CUTOFF
        )


    U_ids <- unique(
        U$ID
    )


    N_ids <- unique(
        N$ID
    )


    common <- intersect(
        U_ids,
        N_ids
    )


    U_only <- setdiff(
        U_ids,
        N_ids
    )


    N_only <- setdiff(
        N_ids,
        U_ids
    )


    all_ids <- unique(
        c(
            U_ids,
            N_ids
        )
    )


    cat(
        "\nUniform significant: ",
        length(U_ids),
        "\n"
    )


    cat(
        "Non-uniform significant: ",
        length(N_ids),
        "\n"
    )


    cat(
        "Common: ",
        length(common),
        "\n"
    )


    cat(
        "Uniform-specific: ",
        length(U_only),
        "\n"
    )


    cat(
        "Non-uniform-specific: ",
        length(N_only),
        "\n"
    )


    # --------------------------------------------------------
    # Prepare Uniform table
    # --------------------------------------------------------

    U_desc <- U %>%

        select(
            ID,
            Description,
            pvalue,
            p.adjust,
            FDR_BH
        ) %>%

        distinct(
            ID,
            .keep_all = TRUE
        ) %>%

        rename(

            Term_ID =
                ID,

            Uniform_Description =
                Description,

            Uniform_pvalue =
                pvalue,

            Uniform_p_adjust =
                p.adjust,

            Uniform_FDR =
                FDR_BH

        )


    # --------------------------------------------------------
    # Prepare Non-uniform table
    # --------------------------------------------------------

    N_desc <- N %>%

        select(
            ID,
            Description,
            pvalue,
            p.adjust,
            FDR_BH
        ) %>%

        distinct(
            ID,
            .keep_all = TRUE
        ) %>%

        rename(

            Term_ID =
                ID,

            NonUniform_Description =
                Description,

            NonUniform_pvalue =
                pvalue,

            NonUniform_p_adjust =
                p.adjust,

            NonUniform_FDR =
                FDR_BH

        )


    # --------------------------------------------------------
    # Master comparison
    # --------------------------------------------------------

    comparison <- data.frame(

        Term_ID =
            all_ids,

        stringsAsFactors =
            FALSE

    )


    comparison <- comparison %>%

        left_join(
            U_desc,
            by = "Term_ID"
        ) %>%

        left_join(
            N_desc,
            by = "Term_ID"
        )


    comparison$Classification <- ifelse(

        comparison$Term_ID %in%
            common,

        "Common",

        ifelse(

            comparison$Term_ID %in%
                U_only,

            "Uniform_specific",

            "NonUniform_specific"

        )

    )


    comparison <- comparison %>%

        select(

            Term_ID,

            Uniform_Description,
            Uniform_pvalue,
            Uniform_p_adjust,
            Uniform_FDR,

            NonUniform_Description,
            NonUniform_pvalue,
            NonUniform_p_adjust,
            NonUniform_FDR,

            Classification

        )


    # --------------------------------------------------------
    # Save comparison
    # --------------------------------------------------------

    comparison_file <- file.path(

        OUTPUT_DIR,

        "Uniform_vs_NonUniform",

        database,

        paste0(
            database,
            "_Uniform_vs_NonUniform.csv"
        )

    )


    write.csv(
        comparison,
        comparison_file,
        row.names = FALSE
    )


    # --------------------------------------------------------
    # Save individual groups
    # --------------------------------------------------------

    write.csv(

        U,

        file.path(

            OUTPUT_DIR,

            "Uniform_vs_NonUniform",

            database,

            "Uniform_SIGNIFICANT.csv"

        ),

        row.names = FALSE

    )


    write.csv(

        N,

        file.path(

            OUTPUT_DIR,

            "Uniform_vs_NonUniform",

            database,

            "NonUniform_SIGNIFICANT.csv"

        ),

        row.names = FALSE

    )


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary <- data.frame(

        Database =
            database,

        Uniform_Significant =
            length(U_ids),

        NonUniform_Significant =
            length(N_ids),

        Common =
            length(common),

        Uniform_Specific =
            length(U_only),

        NonUniform_Specific =
            length(N_only),

        stringsAsFactors =
            FALSE

    )


    write.csv(

        summary,

        file.path(

            OUTPUT_DIR,

            "Uniform_vs_NonUniform",

            database,

            "Summary.csv"

        ),

        row.names = FALSE

    )


    # ========================================================
    # ELLIPSE-SHAPED VENN DIAGRAM
    # ========================================================

    venn_file <- file.path(

        OUTPUT_DIR,

        "Plots",

        paste0(
            database,
            "_Uniform_vs_NonUniform_Ellipse_Venn.png"
        )

    )


    png(

        venn_file,

        width = 3000,

        height = 2400,

        res = 300

    )


    if (
        length(U_ids) > 0 ||
        length(N_ids) > 0
    ) {


        venn_sets <- list(

            "UUni\n(Uniform)" =
                U_ids,

            "UVar\n(Non-uniform)" =
                N_ids

        )


        venn(

            venn_sets,

            # ------------------------------------------------
            # FORCE ELLIPSE SHAPE
            # ------------------------------------------------

            ellipse = TRUE,

            # ------------------------------------------------
            # Show intersection numbers
            # ------------------------------------------------

            ilabels = TRUE,

            # ------------------------------------------------
            # Pink / Grey
            #
            # UUni = Pink
            # UVar = Grey
            # Common = Light grey
            # ------------------------------------------------

            zcolor = c(
                "#F4A6C1",
                "#BDBDBD",
                "#D9D9D9"
            ),

            opacity = 0.65,

            # ------------------------------------------------
            # Borders
            # ------------------------------------------------

            borders = TRUE,

            col = c(
                "#D96B91",
                "#7A7A7A"
            ),

            lwd = 3,

            # ------------------------------------------------
            # No surrounding box
            # ------------------------------------------------

            box = FALSE,

            # ------------------------------------------------
            # Intersection number size
            # ------------------------------------------------

            ilcs = 1.6,

            # ------------------------------------------------
            # Set name size
            # ------------------------------------------------

            sncs = 1.4,

            # ------------------------------------------------
            # Plot size
            # ------------------------------------------------

            plotsize = 18

        )


        # ----------------------------------------------------
        # Bold title
        # ----------------------------------------------------

        grid::grid.text(

            paste0(
                "MIsoTEx ",
                database,
                "\nUbiquitous : Uniform (UUni) vs Non-uniform (UVar)",
                "\nFDR < 0.05"
            ),

            x = 0.5,

            y = 0.96,

            gp = grid::gpar(

                fontsize = 18,

                fontface = "bold"

            )

        )


    } else {


        grid::grid.newpage()

        grid::grid.text(

            "No significant enriched functions",

            gp = grid::gpar(

                fontsize = 20,

                fontface = "bold"

            )

        )

    }


    dev.off()


    cat(
        "\nEllipse-shaped Venn created:\n",
        venn_file,
        "\n"
    )


    # --------------------------------------------------------
    # RETURN
    # --------------------------------------------------------

    list(

        comparison =
            comparison,

        summary =
            summary,

        uniform_ids =
            U_ids,

        nonuniform_ids =
            N_ids,

        common =
            common,

        uniform_specific =
            U_only,

        nonuniform_specific =
            N_only

    )

}


# ============================================================
# 24. RUN ALL COMPARISONS
# ============================================================

COMPARE_GO_BP <- compare_functions(
    UNIFORM_GO_BP,
    NONUNIFORM_GO_BP,
    "GO_BP"
)


COMPARE_GO_MF <- compare_functions(
    UNIFORM_GO_MF,
    NONUNIFORM_GO_MF,
    "GO_MF"
)


COMPARE_GO_CC <- compare_functions(
    UNIFORM_GO_CC,
    NONUNIFORM_GO_CC,
    "GO_CC"
)


COMPARE_KEGG <- compare_functions(
    UNIFORM_KEGG,
    NONUNIFORM_KEGG,
    "KEGG"
)


COMPARE_REACTOME <- compare_functions(
    UNIFORM_REACTOME,
    NONUNIFORM_REACTOME,
    "Reactome"
)


# ============================================================
# 25. OVERALL SUMMARY TABLE
# ============================================================

overall_summary <- bind_rows(

    COMPARE_GO_BP$summary,

    COMPARE_GO_MF$summary,

    COMPARE_GO_CC$summary,

    COMPARE_KEGG$summary,

    COMPARE_REACTOME$summary

)


write.csv(

    overall_summary,

    file.path(

        OUTPUT_DIR,

        "Uniform_vs_NonUniform",

        "FINAL_Uniform_vs_NonUniform_Summary.csv"

    ),

    row.names = FALSE

)


# ============================================================
# 26. CREATE COMBINED SUMMARY IMAGE
# ============================================================

plot_summary <- rbind(

    data.frame(

        Database =
            overall_summary$Database,

        Classification =
            "Uniform_Specific",

        Number =
            overall_summary$Uniform_Specific,

        stringsAsFactors =
            FALSE

    ),

    data.frame(

        Database =
            overall_summary$Database,

        Classification =
            "Common",

        Number =
            overall_summary$Common,

        stringsAsFactors =
            FALSE

    ),

    data.frame(

        Database =
            overall_summary$Database,

        Classification =
            "NonUniform_Specific",

        Number =
            overall_summary$NonUniform_Specific,

        stringsAsFactors =
            FALSE

    )

)


plot_summary$Classification <-

    factor(

        plot_summary$Classification,

        levels = c(

            "Uniform_Specific",

            "Common",

            "NonUniform_Specific"

        )

    )


p_summary <- ggplot(

    plot_summary,

    aes(

        x =
            Database,

        y =
            Number,

        fill =
            Classification

    )

) +

    geom_col(
        position = "dodge"
    ) +

    geom_text(

        aes(
            label = Number
        ),

        position =
            position_dodge(
                width = 0.9
            ),

        vjust = -0.25,

        size = 3.5,

        fontface = "bold"

    ) +

    scale_fill_manual(

        values = c(

            "Uniform_Specific" =
                "#F4A6C1",

            "Common" =
                "#D9D9D9",

            "NonUniform_Specific" =
                "#BDBDBD"

        ),

        labels = c(

            "Uniform-specific",

            "Common",

            "Non-uniform-specific"

        )

    ) +

    labs(

        title =
            "MIsoTEx Functional Comparison",

        subtitle =
            "Pooled Uniform vs Pooled Non-uniform",

        x =
            NULL,

        y =
            "Number of Significant Functions",

        fill =
            NULL

    ) +

    theme_minimal(

        base_size = 13

    ) +

    theme(

        axis.text.x =
            element_text(
                angle = 30,
                hjust = 1,
                face = "bold"
            ),

        axis.text.y =
            element_text(
                face = "bold"
            ),

        axis.title.y =
            element_text(
                face = "bold"
            ),

        plot.title =
            element_text(
                face = "bold",
                size = 17
            ),

        plot.subtitle =
            element_text(
                face = "bold"
            ),

        legend.text =
            element_text(
                face = "bold"
            ),

        legend.position =
            "bottom"

    )


ggsave(

    file.path(

        OUTPUT_DIR,

        "Plots",

        "FINAL_Uniform_vs_NonUniform_Summary.png"

    ),

    p_summary,

    width = 12,

    height = 7,

    dpi = 300

)


# ============================================================
# 27. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("MIsoTEx FUNCTIONAL ANALYSIS COMPLETE\n")
cat("============================================================\n\n")


cat(
    "IMPORTANT ANALYSIS DESIGN:\n\n"
)


cat(
    "Uniform = ALL expression types pooled\n"
)


cat(
    "Non-uniform = ALL expression types pooled\n\n"
)


cat(
    "Very High + High + Moderate + Low + Very Low\n"
)


cat(
    "were NOT compared separately.\n\n"
)


cat(
    "No >50 transcript filtering was applied.\n"
)


cat(
    "Multiple-testing correction: Benjamini-Hochberg\n"
)


cat(
    "Significance cutoff: FDR < ",
    FDR_CUTOFF,
    "\n\n"
)


cat(
    "Venn diagram: ELLIPSE-SHAPED\n"
)


cat(
    "Uniform / UUni: PINK\n"
)


cat(
    "Non-uniform / UVar: GREY\n\n"
)


cat(
    "Uniform transcripts: ",
    nrow(uniform_df),
    "\n"
)


cat(
    "Non-uniform transcripts: ",
    nrow(nonuniform_df),
    "\n\n"
)


cat(
    "Final output directory:\n",
    OUTPUT_DIR,
    "\n\n"
)


cat(
    "Main comparison summary:\n",
    file.path(
        OUTPUT_DIR,
        "Uniform_vs_NonUniform",
        "FINAL_Uniform_vs_NonUniform_Summary.csv"
    ),
    "\n\n"
)


cat(
    "Final summary image:\n",
    file.path(
        OUTPUT_DIR,
        "Plots",
        "FINAL_Uniform_vs_NonUniform_Summary.png"
    ),
    "\n\n"
)


cat(
    "Ellipse Venn diagrams:\n",
    file.path(
        OUTPUT_DIR,
        "Plots"
    ),
    "\n\n"
)


cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")
