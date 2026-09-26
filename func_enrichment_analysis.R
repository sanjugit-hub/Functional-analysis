# ------------------------------------------------------------
# Purpose:
# Perform functional enrichment analysis of MIsoTEx ubiquitous
# transcripts classified into five expression levels.
#
# Analyses:
# GO-BP, GO-MF, GO-CC, KEGG and Reactome enrichment,
# followed by function overlap, Venn and colored UpSet plots.
#
# Author: SBM
# ------------------------------------------------------------


# ============================================================
# 1. USER SETTINGS
# ============================================================

#INPUT_FILE <- "Varying_quantity_of_ubiquitously_present_list_Expression_profile_Score.csv"
INPUT_FILE <- "testis_iso-mRNA_Predominantly_present_in_tissues_Expression_profile_Score.csv"
OUTPUT_DIR <- "MIsoTEx_Testis_Functional_Analysis"

ADJ_P_CUTOFF <- 0.05

TOP_N <- 20

MAX_INTERSECTIONS <- 30


# ============================================================
# 2. EXPRESSION LEVEL ORDER
# ============================================================

expression_levels <- c(
    "Very High",
    "High",
    "Moderate",
    "Low",
    "Very Low"
)


# ============================================================
# 3. COLOUR PALETTE
# ============================================================

# Very High = dark red
# High      = pink/red
# Moderate  = orange
# Low       = green
# Very Low  = purple

expression_colors <- c(

    "Very High" = "#990000",

    "High" = "#F08080",

    "Moderate" = "#FFA500",

    "Low" = "#228B22",

    "Very Low" = "#6A5ACD"

)

intersection_color <- "#4682B4"


# ============================================================
# 4. REQUIRED PACKAGES
# ============================================================

required_packages <- c(

    "clusterProfiler",

    "enrichplot",

    "KEGGREST",

    "ReactomePA",

    "org.Hs.eg.db",

    "AnnotationDbi",

    "VennDiagram",

    "ggplot2",

    "dplyr",

    "grid"

)


# ============================================================
# 5. CHECK PACKAGES
# ============================================================

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
        "This script does NOT install packages automatically.\n"

    )

}


cat("All required packages are installed.\n\n")


# ============================================================
# 6. LOAD PACKAGES
# ============================================================

suppressPackageStartupMessages({

    library(clusterProfiler)

    library(enrichplot)

    library(KEGGREST)

    library(ReactomePA)

    library(org.Hs.eg.db)

    library(AnnotationDbi)

    library(VennDiagram)

    library(ggplot2)

    library(dplyr)

    library(grid)

})


# ============================================================
# 7. PACKAGE VERSIONS
# ============================================================

cat("PACKAGE VERSIONS\n")
cat("----------------\n\n")


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
# 8. CREATE OUTPUT DIRECTORIES
# ============================================================

output_dirs <- c(

    "",

    "Individual_GO",

    "Individual_KEGG",

    "Individual_Reactome",

    "GO_BP",

    "GO_MF",

    "GO_CC",

    "KEGG",

    "Reactome",

    "Comparison",

    "Plots",

    "Plots/GO_BP",

    "Plots/GO_MF",

    "Plots/GO_CC",

    "Plots/KEGG",

    "Plots/Reactome",

    "Transcript_Mapping"

)


for (d in output_dirs) {

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
# 9. CHECK INPUT FILE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CHECKING INPUT FILE\n")
cat("============================================================\n\n")


if (!file.exists(INPUT_FILE)) {

    stop(

        "\nInput file not found:\n",

        INPUT_FILE,

        "\n\nCurrent working directory:\n",

        getwd(),

        "\n"

    )

}


cat(
    "Input file: ",
    INPUT_FILE,
    "\n"
)


# ============================================================
# 10. READ INPUT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("READING INPUT FILE\n")
cat("============================================================\n\n")


df <- read.csv(

    INPUT_FILE,

    header = TRUE,

    stringsAsFactors = FALSE,

    check.names = FALSE

)


cat(
    "Total rows: ",
    nrow(df),
    "\n\n"
)


cat("Input columns:\n")

print(
    colnames(df)
)


# ============================================================
# 11. STANDARDIZE COLUMN NAMES
# ============================================================

name_map <- c(

    "Gene_stable_ID" = "gene_stable_id",

    "Gene_Stable_ID" = "gene_stable_id",

    "gene_stable_ID" = "gene_stable_id",

    "Gene_name" = "gene_name",

    "Gene_Name" = "gene_name",

    "Expression_Type" = "Expression_type"

)


for (old_name in names(name_map)) {

    if (old_name %in% colnames(df)) {

        colnames(df)[

            colnames(df) == old_name

        ] <- name_map[old_name]

    }

}


cat("\nStandardized columns:\n")

print(
    colnames(df)
)


# ============================================================
# 12. REQUIRED COLUMN CHECK
# ============================================================

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

        "\nMissing required columns:\n",

        paste(

            missing_columns,

            collapse = ", "

        ),

        "\n"

    )

}


# ============================================================
# 13. CLEAN IDs
# ============================================================

df$transcript_id <- sub(

    "\\..*$",

    "",

    as.character(
        df$transcript_id
    )

)


df$gene_stable_id <- sub(

    "\\..*$",

    "",

    as.character(
        df$gene_stable_id
    )

)


df$Expression_type <- trimws(

    as.character(
        df$Expression_type
    )

)


df <- df %>%

    filter(

        !is.na(
            transcript_id
        ),

        transcript_id != "",

        !is.na(
            gene_stable_id
        ),

        gene_stable_id != "",

        !is.na(
            Expression_type
        ),

        Expression_type != ""

    )


cat(

    "\nRows after cleaning: ",

    nrow(df),

    "\n"

)


# ============================================================
# 14. EXPRESSION LEVEL CHECK
# ============================================================

cat("\n")
cat("============================================================\n")
cat("EXPRESSION LEVEL DISTRIBUTION\n")
cat("============================================================\n\n")


print(
    table(
        df$Expression_type
    )
)


unexpected_levels <- setdiff(

    unique(
        df$Expression_type
    ),

    expression_levels

)


if (length(unexpected_levels) > 0) {

    warning(

        "\nUnexpected expression levels found:\n",

        paste(

            unexpected_levels,

            collapse = ", "

        )

    )

}


# ============================================================
# 15. EXPRESSION SUMMARY
# ============================================================

expression_summary <- data.frame(

    Expression_Type =
        expression_levels,

    Transcript_Count =
        sapply(

            expression_levels,

            function(x) {

                sum(
                    df$Expression_type == x
                )

            }

        ),

    Unique_Gene_Count =
        sapply(

            expression_levels,

            function(x) {

                length(

                    unique(

                        df$gene_stable_id[

                            df$Expression_type == x

                        ]

                    )

                )

            }

        ),

    stringsAsFactors = FALSE

)


write.csv(

    expression_summary,

    file.path(

        OUTPUT_DIR,

        "Expression_Level_Summary.csv"

    ),

    row.names = FALSE

)


cat("\nExpression summary:\n\n")

print(
    expression_summary
)


# ============================================================
# 16. ENSEMBL → ENTREZ MAPPING
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ENSEMBL TO ENTREZ MAPPING\n")
cat("============================================================\n\n")


all_genes <- unique(

    df$gene_stable_id

)


cat(

    "Unique Ensembl genes: ",

    length(all_genes),

    "\n"

)


gene_conversion <- tryCatch(

    {

        clusterProfiler::bitr(

            all_genes,

            fromType = "ENSEMBL",

            toType = "ENTREZID",

            OrgDb = org.Hs.eg.db

        )

    },

    error = function(e) {

        stop(

            "\nEnsembl to Entrez mapping failed:\n",

            e$message,

            "\n"

        )

    }

)


gene_conversion <- gene_conversion %>%

    distinct(

        ENSEMBL,

        ENTREZID,

        .keep_all = TRUE

    )


mapped_genes <- length(

    unique(
        gene_conversion$ENSEMBL
    )

)


mapping_percentage <- round(

    100 *

    mapped_genes /

    length(all_genes),

    2

)


cat(

    "Mapped genes: ",

    mapped_genes,

    "\n"

)


cat(

    "Mapping percentage: ",

    mapping_percentage,

    "%\n"

)


write.csv(

    gene_conversion,

    file.path(

        OUTPUT_DIR,

        "Transcript_Mapping",

        "Ensembl_to_Entrez.csv"

    ),

    row.names = FALSE

)


# ============================================================
# 17. TRANSCRIPT → GENE → ENTREZ MAPPING
# ============================================================

transcript_mapping <- df %>%

    select(

        transcript_id,

        gene_stable_id,

        gene_name,

        Expression_type

    ) %>%

    left_join(

        gene_conversion,

        by = c(

            "gene_stable_id" = "ENSEMBL"

        )

    )


write.csv(

    transcript_mapping,

    file.path(

        OUTPUT_DIR,

        "Transcript_Mapping",

        "Transcript_Gene_Entrez_Mapping.csv"

    ),

    row.names = FALSE

)


# ============================================================
# 18. COMPLETE HUMAN BACKGROUND
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PREPARING HUMAN BACKGROUND\n")
cat("============================================================\n\n")


background_entrez <- AnnotationDbi::keys(

    org.Hs.eg.db,

    keytype = "ENTREZID"

)


background_entrez <- unique(

    background_entrez[

        !is.na(
            background_entrez
        )

    ]

)


cat(

    "Background Entrez genes: ",

    length(background_entrez),

    "\n"

)


# ============================================================
# 19. TRANSCRIPT CONTRIBUTOR FUNCTION
# ============================================================

get_transcript_contributors <- function(

    gene_string,

    level_df

) {


    if (

        is.na(gene_string) ||

        gene_string == ""

    ) {

        return("")

    }


    entrez_ids <- unique(

        unlist(

            strsplit(

                gene_string,

                "/"

            )

        )

    )


    contributors <- level_df %>%

        filter(

            ENTREZID %in% entrez_ids

        )


    paste(

        unique(

            contributors$transcript_id

        ),

        collapse = ";"

    )

}


# ============================================================
# 20. ADD STATISTICAL COLUMNS
# ============================================================

add_statistics <- function(

    result,

    level

) {


    result$Expression_Type <- level


    # BH adjusted p-value = FDR

    result$FDR_BH <- result$p.adjust


    result$Significant <- ifelse(

        !is.na(
            result$FDR_BH
        ) &

        result$FDR_BH < ADJ_P_CUTOFF,

        "Yes",

        "No"

    )


    result

}


# ============================================================
# 21. GO ENRICHMENT FUNCTION
# ============================================================

run_GO <- function(

    level,

    ontology

) {


    cat("\n")

    cat(

        "Running GO-",

        ontology,

        " : ",

        level,

        "\n",

        sep = ""

    )


    level_df <- transcript_mapping %>%

        filter(

            Expression_type == level

        )


    transcript_ids <- unique(

        level_df$transcript_id

    )


    genes <- unique(

        level_df$ENTREZID[

            !is.na(
                level_df$ENTREZID
            )

        ]

    )


    cat(

        "Transcripts: ",

        length(transcript_ids),

        "\n"

    )


    cat(

        "Mapped unique genes: ",

        length(genes),

        "\n"

    )


    if (length(genes) == 0) {

        warning(

            "No mapped genes for ",

            level

        )

        return(NULL)

    }


    ego <- tryCatch(

        {

            clusterProfiler::enrichGO(

                gene = genes,

                universe = background_entrez,

                OrgDb = org.Hs.eg.db,

                keyType = "ENTREZID",

                ont = ontology,

                pAdjustMethod = "BH",

                pvalueCutoff = 1,

                qvalueCutoff = 1,

                readable = FALSE

            )

        },

        error = function(e) {

            warning(

                "GO failed for ",

                level,

                " / ",

                ontology,

                ": ",

                e$message

            )

            NULL

        }

    )


    if (is.null(ego)) {

        return(NULL)

    }


    result <- as.data.frame(
        ego
    )


    if (nrow(result) == 0) {

        return(NULL)

    }


    result <- add_statistics(

        result,

        level

    )


    result$Ontology <- ontology


    result$Transcript_IDs <- sapply(

        result$geneID,

        get_transcript_contributors,

        level_df = level_df

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


    result

}


# ============================================================
# 22. RUN GO
# ============================================================

cat("\n")
cat("============================================================\n")
cat("RUNNING GO ENRICHMENT\n")
cat("============================================================\n")


GO_RESULTS <- list()


for (ontology in c(

    "BP",

    "MF",

    "CC"

)) {


    for (level in expression_levels) {


        key <- paste(

            ontology,

            level,

            sep = "_"

        )


        GO_RESULTS[[key]] <- run_GO(

            level,

            ontology

        )

    }

}


# ============================================================
# 23. SAVE GO RESULTS
# ============================================================

for (ontology in c(

    "BP",

    "MF",

    "CC"

)) {


    combined <- list()


    for (level in expression_levels) {


        key <- paste(

            ontology,

            level,

            sep = "_"

        )


        result <- GO_RESULTS[[key]]


        if (is.null(result)) {

            next

        }


        combined[[level]] <- result


        safe_level <- gsub(

            " ",

            "_",

            level

        )


        write.csv(

            result,

            file.path(

                OUTPUT_DIR,

                "Individual_GO",

                paste0(

                    "GO_",

                    ontology,

                    "_",

                    safe_level,

                    "_ALL.csv"

                )

            ),

            row.names = FALSE

        )


        write.csv(

            result %>%

                filter(

                    FDR_BH < ADJ_P_CUTOFF

                ),

            file.path(

                OUTPUT_DIR,

                "Individual_GO",

                paste0(

                    "GO_",

                    ontology,

                    "_",

                    safe_level,

                    "_SIGNIFICANT.csv"

                )

            ),

            row.names = FALSE

        )

    }


    if (length(combined) > 0) {


        combined_df <- bind_rows(

            combined

        )


        write.csv(

            combined_df,

            file.path(

                OUTPUT_DIR,

                paste0(

                    "GO_",

                    ontology

                ),

                "ALL_LEVELS.csv"

            ),

            row.names = FALSE

        )


        write.csv(

            combined_df %>%

                filter(

                    FDR_BH < ADJ_P_CUTOFF

                ),

            file.path(

                OUTPUT_DIR,

                paste0(

                    "GO_",

                    ontology

                ),

                "SIGNIFICANT_ALL_LEVELS.csv"

            ),

            row.names = FALSE

        )

    }

}


# ============================================================
# 24. FUNCTION PRESENCE / ABSENCE MATRIX
# ============================================================

create_function_matrix <- function(

    result_list,

    database_name

) {


    term_sets <- list()


    for (level in expression_levels) {


        result <- result_list[[level]]


        if (is.null(result)) {

            term_sets[[level]] <-
                character(0)

        } else {

            term_sets[[level]] <- unique(

                result$ID[

                    result$FDR_BH < ADJ_P_CUTOFF

                ]

            )

        }

    }


    all_terms <- sort(

        unique(

            unlist(
                term_sets
            )

        )

    )


    if (length(all_terms) == 0) {

        cat(

            "\nNo significant functions for ",

            database_name,

            "\n"

        )

        return(NULL)

    }


    matrix_df <- data.frame(

        Term_ID = all_terms,

        stringsAsFactors = FALSE

    )


    for (level in expression_levels) {

        matrix_df[[level]] <- as.integer(

            all_terms %in%

            term_sets[[level]]

        )

    }


    matrix_df$Number_of_Levels <- rowSums(

        matrix_df[

            ,

            expression_levels

        ]

    )


    # --------------------------------------------------------
    # Add description
    # --------------------------------------------------------

    descriptions <- list()


    for (level in expression_levels) {


        result <- result_list[[level]]


        if (!is.null(result)) {

            descriptions[[level]] <-

                result %>%

                select(

                    ID,

                    Description

                ) %>%

                distinct()

        }

    }


    descriptions <- bind_rows(

        descriptions

    ) %>%

        distinct(

            ID,

            .keep_all = TRUE

        )


    matrix_df <- matrix_df %>%

        left_join(

            descriptions,

            by = c(

                "Term_ID" = "ID"

            )

        ) %>%

        select(

            Term_ID,

            Description,

            everything()

        )


    write.csv(

        matrix_df,

        file.path(

            OUTPUT_DIR,

            "Comparison",

            paste0(

                database_name,

                "_Function_Presence_Absence_Matrix.csv"

            )

        ),

        row.names = FALSE

    )


    matrix_df

}


GO_BP_MATRIX <- create_function_matrix(

    setNames(

        lapply(

            expression_levels,

            function(x) {

                GO_RESULTS[[

                    paste(

                        "BP",

                        x,

                        sep = "_"

                    )

                ]]

            }

        ),

        expression_levels

    ),

    "GO_BP"

)


GO_MF_MATRIX <- create_function_matrix(

    setNames(

        lapply(

            expression_levels,

            function(x) {

                GO_RESULTS[[

                    paste(

                        "MF",

                        x,

                        sep = "_"

                    )

                ]]

            }

        ),

        expression_levels

    ),

    "GO_MF"

)


GO_CC_MATRIX <- create_function_matrix(

    setNames(

        lapply(

            expression_levels,

            function(x) {

                GO_RESULTS[[

                    paste(

                        "CC",

                        x,

                        sep = "_"

                    )

                ]]

            }

        ),

        expression_levels

    ),

    "GO_CC"

)


# ============================================================
# 25. OVERLAP SUMMARY
# ============================================================

create_overlap_summary <- function(

    matrix_df,

    database_name

) {


    if (is.null(matrix_df)) {

        return(NULL)

    }


    n <- matrix_df$Number_of_Levels


    summary_df <- data.frame(

        Category = c(

            "Common to all 5 levels",

            "Common to exactly 4 levels",

            "Common to exactly 3 levels",

            "Common to exactly 2 levels",

            "Unique to one level",

            "Total unique enriched functions"

        ),

        Number_of_Functions = c(

            sum(n == 5),

            sum(n == 4),

            sum(n == 3),

            sum(n == 2),

            sum(n == 1),

            nrow(matrix_df)

        )

    )


    write.csv(

        summary_df,

        file.path(

            OUTPUT_DIR,

            "Comparison",

            paste0(

                database_name,

                "_Overlap_Summary.csv"

            )

        ),

        row.names = FALSE

    )


    # --------------------------------------------------------
    # Exact combinations
    # --------------------------------------------------------

    exact_combination <- apply(

        matrix_df[

            ,

            expression_levels

        ],

        1,

        function(x) {


            present <- expression_levels[

                x == 1

            ]


            if (length(present) == 0) {

                return("None")

            }


            paste(

                present,

                collapse = " + "

            )

        }

    )


    exact_df <- data.frame(

        Exact_Combination =
            exact_combination,

        stringsAsFactors = FALSE

    ) %>%

        count(

            Exact_Combination,

            name = "Number_of_Functions"

        ) %>%

        arrange(

            desc(
                Number_of_Functions
            )

        )


    write.csv(

        exact_df,

        file.path(

            OUTPUT_DIR,

            "Comparison",

            paste0(

                database_name,

                "_Exact_Intersection_Combinations.csv"

            )

        ),

        row.names = FALSE

    )


    list(

        summary = summary_df,

        exact = exact_df

    )

}


BP_SUMMARY <- create_overlap_summary(

    GO_BP_MATRIX,

    "GO_BP"

)


MF_SUMMARY <- create_overlap_summary(

    GO_MF_MATRIX,

    "GO_MF"

)


CC_SUMMARY <- create_overlap_summary(

    GO_CC_MATRIX,

    "GO_CC"

)


# ============================================================
# 26. COLOURED UPSET PLOT
# ============================================================

create_custom_upset <- function(

    matrix_df,

    database_name

) {


    cat("\n")
    cat("============================================================\n")
    cat("COLOURED UPSET PLOT: ", database_name, "\n")
    cat("============================================================\n")


    if (is.null(matrix_df)) {

        cat(
            "No matrix available.\n"
        )

        return(NULL)

    }


    if (nrow(matrix_df) == 0) {

        cat(
            "Matrix contains zero functions.\n"
        )

        return(NULL)

    }


    # --------------------------------------------------------
    # Set sizes
    # --------------------------------------------------------

    set_sizes <- sapply(

        expression_levels,

        function(x) {

            sum(

                matrix_df[[x]] == 1,

                na.rm = TRUE

            )

        }

    )


    cat(
        "\nSignificant functions per expression level:\n\n"
    )


    print(
        set_sizes
    )


    # --------------------------------------------------------
    # Keep non-empty levels
    # --------------------------------------------------------

    active_levels <- expression_levels[

        set_sizes > 0

    ]


    if (length(active_levels) < 2) {

        cat(
            "\nFewer than two active sets.\n"
        )

        return(NULL)

    }


    # --------------------------------------------------------
    # Binary matrix
    # --------------------------------------------------------

    m <- matrix_df %>%

        select(

            all_of(
                active_levels
            )

        )


    m[] <- lapply(

        m,

        function(x) {

            as.integer(x)

        }

    )


    # --------------------------------------------------------
    # Remove functions with no membership
    # --------------------------------------------------------

    m <- m[

        rowSums(m) > 0,

        ,

        drop = FALSE

    ]


    # --------------------------------------------------------
    # Intersection pattern
    # --------------------------------------------------------

    pattern <- apply(

        m,

        1,

        function(x) {

            paste0(

                x,

                collapse = ""

            )

        }

    )


    pattern_count <- sort(

        table(pattern),

        decreasing = TRUE

    )


    intersection_df <- data.frame(

        Pattern =
            names(pattern_count),

        Intersection_Size =
            as.integer(pattern_count),

        stringsAsFactors = FALSE

    )


    # --------------------------------------------------------
    # Maximum displayed intersections
    # --------------------------------------------------------

    if (

        nrow(intersection_df) >

        MAX_INTERSECTIONS

    ) {

        intersection_df <-

            intersection_df[

                1:MAX_INTERSECTIONS,

                ,

                drop = FALSE

            ]

    }


    intersection_df$Intersection <-

        seq_len(

            nrow(intersection_df)

        )


    # --------------------------------------------------------
    # Save intersection table
    # --------------------------------------------------------

    safe_name <- gsub(

        "[^A-Za-z0-9]+",

        "_",

        database_name

    )


    write.csv(

        intersection_df,

        file.path(

            OUTPUT_DIR,

            "Comparison",

            paste0(

                safe_name,

                "_UpSet_Intersection_Summary.csv"

            )

        ),

        row.names = FALSE

    )


    # ========================================================
    # DOT MATRIX DATA
    # ========================================================

    all_dot_df <- expand.grid(

        Intersection =

            intersection_df$Intersection,

        Expression_Level =

            active_levels,

        stringsAsFactors = FALSE

    )


    all_dot_df$Pattern <- mapply(

        function(intersection_number) {

            intersection_df$Pattern[

                intersection_df$Intersection ==

                    intersection_number

            ]

        },

        all_dot_df$Intersection

    )


    all_dot_df$Position <- match(

        all_dot_df$Expression_Level,

        active_levels

    )


    all_dot_df$Present <- mapply(

        function(

            pattern_value,

            position

        ) {

            as.integer(

                substr(

                    pattern_value,

                    position,

                    position

                )

            )

        },

        all_dot_df$Pattern,

        all_dot_df$Position

    )


    # --------------------------------------------------------
    # IMPORTANT:
    # Reverse factor levels so plot appears:
    #
    # Very High
    # High
    # Moderate
    # Low
    # Very Low
    #
    # from top to bottom.
    # --------------------------------------------------------

    all_dot_df$Expression_Level <-

        factor(

            all_dot_df$Expression_Level,

            levels =
                rev(active_levels)

        )


    # ========================================================
    # CONNECTING LINE DATA
    # ========================================================

    line_df <- data.frame()


    for (i in seq_len(

        nrow(intersection_df)

    )) {


        pattern_value <-

            intersection_df$Pattern[i]


        values <- as.integer(

            strsplit(

                pattern_value,

                split = ""

            )[[1]]

        )


        present_positions <- which(

            values == 1

        )


        if (

            length(
                present_positions
            ) > 1

        ) {


            line_df <- rbind(

                line_df,

                data.frame(

                    Intersection =

                        intersection_df$

                        Intersection[i],

                    ymin =

                        active_levels[

                            min(
                                present_positions
                            )

                        ],

                    ymax =

                        active_levels[

                            max(
                                present_positions
                            )

                        ]

                )

            )

        }

    }


    if (nrow(line_df) > 0) {

        line_df$ymin <-

            factor(

                line_df$ymin,

                levels =
                    rev(active_levels)

            )


        line_df$ymax <-

            factor(

                line_df$ymax,

                levels =
                    rev(active_levels)

            )

    }


    # ========================================================
    # SET SIZE DATA
    # ========================================================

    set_size_df <- data.frame(

        Expression_Level =

            active_levels,

        Set_Size =

            as.integer(

                set_sizes[
                    active_levels
                ]

            ),

        stringsAsFactors = FALSE

    )


    set_size_df$Expression_Level <-

        factor(

            set_size_df$Expression_Level,

            levels =
                rev(active_levels)

        )


    # ========================================================
    # TOP INTERSECTION BAR
    # ========================================================

    p_bar <- ggplot(

        intersection_df,

        aes(

            x = Intersection,

            y = Intersection_Size

        )

    ) +

        geom_col(

            width = 0.70,

            fill = intersection_color

        ) +

        geom_text(

            aes(

                label =
                    Intersection_Size

            ),

            vjust = -0.35,

            size = 3.2,

            fontface = "bold"

        ) +

        scale_x_continuous(

            breaks =
                intersection_df$Intersection,

            labels =
                intersection_df$Intersection

        ) +

        scale_y_continuous(

            expand = expansion(

                mult = c(
                    0,
                    0.15
                )

            )

        ) +

        labs(

            x = NULL,

            y = "Intersection Size"

        ) +

        theme_minimal(

            base_size = 12

        ) +

        theme(

            panel.grid.major.x =
                element_blank(),

            panel.grid.minor =
                element_blank(),

            axis.text.x =
                element_blank(),

            axis.ticks.x =
                element_blank(),

            axis.title.y =
                element_text(
                    face = "bold"
                )

        )


    # ========================================================
    # COLOURED DOT MATRIX
    # ========================================================

    p_matrix <- ggplot(

        all_dot_df,

        aes(

            x = Intersection,

            y = Expression_Level

        )

    ) +


        # ----------------------------------------------------
        # Empty dots
        # ----------------------------------------------------

        geom_point(

            data = subset(

                all_dot_df,

                Present == 0

            ),

            shape = 21,

            size = 4.2,

            fill = "grey90",

            color = "grey75",

            stroke = 0.8

        ) +


        # ----------------------------------------------------
        # Present dots
        # ----------------------------------------------------

        geom_point(

            data = subset(

                all_dot_df,

                Present == 1

            ),

            aes(

                fill =
                    Expression_Level

            ),

            shape = 21,

            size = 4.8,

            color = "black",

            stroke = 0.8

        ) +


        scale_fill_manual(

            values =
                expression_colors,

            drop = FALSE

        ) +


        scale_x_continuous(

            breaks =
                intersection_df$Intersection,

            labels =
                intersection_df$Intersection

        ) +


        scale_y_discrete(

            limits =
                rev(active_levels)

        ) +


        labs(

            x = "Intersection",

            y = NULL

        ) +


        theme_minimal(

            base_size = 12

        ) +

        theme(

            panel.grid.major.x =
                element_blank(),

            panel.grid.minor.x =
                element_blank(),

            panel.grid.major.y =
                element_blank(),

            axis.text.y =
                element_text(
                    face = "bold"
                ),

            axis.text.x =
                element_text(
                    size = 10
                ),

            axis.title.x =
                element_text(
                    face = "bold"
                ),

            legend.position =
                "none"

        )


    # --------------------------------------------------------
    # Add black connecting lines
    # --------------------------------------------------------

    if (

        nrow(line_df) > 0

    ) {


        p_matrix <- p_matrix +

            geom_segment(

                data = line_df,

                aes(

                    x =
                        Intersection,

                    xend =
                        Intersection,

                    y =
                        ymin,

                    yend =
                        ymax

                ),

                inherit.aes = FALSE,

                color = "black",

                linewidth = 1.1

            )

    }


    # ========================================================
    # COLOURED SET SIZE BAR
    # ========================================================

    p_set <- ggplot(

        set_size_df,

        aes(

            x = Set_Size,

            y = Expression_Level,

            fill = Expression_Level

        )

    ) +

        geom_col(

            width = 0.65

        ) +

        geom_text(

            aes(

                label =
                    Set_Size

            ),

            hjust = -0.25,

            size = 3.5,

            fontface = "bold"

        ) +

        scale_fill_manual(

            values =
                expression_colors,

            drop = FALSE

        ) +

        scale_y_discrete(

            limits =
                rev(active_levels)

        ) +

        scale_x_continuous(

            expand = expansion(

                mult = c(
                    0,
                    0.18
                )

            )

        ) +

        labs(

            x = "Set Size",

            y = NULL

        ) +

        theme_minimal(

            base_size = 12

        ) +

        theme(

            panel.grid.major.y =
                element_blank(),

            panel.grid.minor =
                element_blank(),

            axis.text.y =
                element_text(
                    face = "bold"
                ),

            axis.title.x =
                element_text(
                    face = "bold"
                ),

            legend.position =
                "none"

        )


    # ========================================================
    # SAVE INDIVIDUAL COMPONENTS
    # ========================================================

    ggsave(

        file.path(

            OUTPUT_DIR,

            "Plots",

            paste0(

                safe_name,

                "_UpSet_Intersection_Bars.png"

            )

        ),

        p_bar,

        width = 12,

        height = 4,

        dpi = 300

    )


    ggsave(

        file.path(

            OUTPUT_DIR,

            "Plots",

            paste0(

                safe_name,

                "_UpSet_Matrix.png"

            )

        ),

        p_matrix,

        width = 12,

        height = 5,

        dpi = 300

    )


    ggsave(

        file.path(

            OUTPUT_DIR,

            "Plots",

            paste0(

                safe_name,

                "_UpSet_Set_Size.png"

            )

        ),

        p_set,

        width = 6,

        height = 5,

        dpi = 300

    )


    # ========================================================
    # COMBINED COLOURED UPSET FIGURE
    # ========================================================

    combined_file <- file.path(

        OUTPUT_DIR,

        "Plots",

        paste0(

            safe_name,

            "_UpSet_Plot.png"

        )

    )


    png(

        combined_file,

        width = 3200,

        height = 2400,

        res = 250

    )


    grid.newpage()


    # --------------------------------------------------------
    # Main title
    # --------------------------------------------------------

    grid.text(

        paste(

            "MIsoTEx -",

            database_name

        ),

        x = 0.5,

        y = 0.975,

        gp = gpar(

            fontsize = 23,

            fontface = "bold"

        )

    )


    grid.text(

        paste(

            "Significantly enriched functions |",

            "BH adjusted p <",

            ADJ_P_CUTOFF

        ),

        x = 0.5,

        y = 0.948,

        gp = gpar(

            fontsize = 13

        )

    )


    # --------------------------------------------------------
    # Top bar
    # --------------------------------------------------------

    print(

        p_bar,

        vp = viewport(

            x = 0.60,

            y = 0.77,

            width = 0.72,

            height = 0.30

        )

    )


    # --------------------------------------------------------
    # Matrix
    # --------------------------------------------------------

    print(

        p_matrix,

        vp = viewport(

            x = 0.60,

            y = 0.47,

            width = 0.72,

            height = 0.34

        )

    )


    # --------------------------------------------------------
    # Set size
    # --------------------------------------------------------

    print(

        p_set,

        vp = viewport(

            x = 0.16,

            y = 0.47,

            width = 0.25,

            height = 0.34

        )

    )


    # --------------------------------------------------------
    # Footer
    # --------------------------------------------------------

    grid.text(

        paste(

            "Total significant functions:",

            nrow(m),

            "|",

            "Intersections displayed:",

            nrow(intersection_df)

        ),

        x = 0.5,

        y = 0.075,

        gp = gpar(

            fontsize = 13,

            fontface = "bold"

        )

    )


    dev.off()


    # ========================================================
    # VERIFY OUTPUT
    # ========================================================

    if (

        file.exists(
            combined_file
        ) &&

        file.info(
            combined_file
        )$size > 1000

    ) {


        cat("\n")

        cat(

            "SUCCESS: Colored UpSet plot created:\n",

            combined_file,

            "\n"

        )

    } else {

        cat(

            "\nWARNING: UpSet plot file was not created correctly.\n"

        )

    }


    invisible(

        list(

            matrix = m,

            intersections =
                intersection_df,

            set_sizes =
                set_size_df

        )

    )

}


# ============================================================
# 27. RUN COLOURED UPSET
# ============================================================

GO_BP_UPSET <- create_custom_upset(

    GO_BP_MATRIX,

    "GO_BP"

)


GO_MF_UPSET <- create_custom_upset(

    GO_MF_MATRIX,

    "GO_MF"

)


GO_CC_UPSET <- create_custom_upset(

    GO_CC_MATRIX,

    "GO_CC"

)


# ============================================================
# 28. VENN DIAGRAM
# ============================================================

create_venn_plot <- function(

    matrix_df,

    database_name

) {


    if (is.null(matrix_df)) {

        return(NULL)

    }


    sets <- list()


    for (level in expression_levels) {

        sets[[level]] <-

            matrix_df$Term_ID[

                matrix_df[[level]] == 1

            ]

    }


    sets <- sets[

        lengths(sets) > 0

    ]


    if (length(sets) < 2) {

        return(NULL)

    }


    output_file <- file.path(

        OUTPUT_DIR,

        "Plots",

        paste0(

            database_name,

            "_Venn.png"

        )

    )


    venn_colors <- c(

        "#990000",

        "#F08080",

        "#FFA500",

        "#228B22",

        "#6A5ACD"

    )


    venn_colors <- venn_colors[

        seq_along(sets)

    ]


    png(

        output_file,

        width = 2600,

        height = 2200,

        res = 250

    )


    venn_plot <- VennDiagram::venn.diagram(

        x = sets,

        filename = NULL,

        fill = venn_colors,

        alpha = 0.45,

        cex = 1,

        cat.cex = 0.9,

        cat.fontface = "bold",

        margin = 0.08,

        main = paste(

            "MIsoTEx - ",

            database_name,

            "\nSignificantly Enriched Functions"

        ),

        main.cex = 1.5

    )


    grid.newpage()

    grid.draw(
        venn_plot
    )


    dev.off()


    cat(

        "\nVenn plot saved:\n",

        output_file,

        "\n"

    )

}


create_venn_plot(

    GO_BP_MATRIX,

    "GO_BP"

)


create_venn_plot(

    GO_MF_MATRIX,

    "GO_MF"

)


create_venn_plot(

    GO_CC_MATRIX,

    "GO_CC"

)


# ============================================================
# 29. VENN + SUMMARY IMAGE
# ============================================================

create_venn_summary <- function(

    matrix_df,

    database_name

) {


    if (is.null(matrix_df)) {

        return(NULL)

    }


    sets <- list()


    for (level in expression_levels) {

        sets[[level]] <-

            matrix_df$Term_ID[

                matrix_df[[level]] == 1

            ]

    }


    sets <- sets[

        lengths(sets) > 0

    ]


    if (length(sets) < 2) {

        return(NULL)

    }


    venn_colors <- c(

        "#990000",

        "#F08080",

        "#FFA500",

        "#228B22",

        "#6A5ACD"

    )


    venn_colors <- venn_colors[

        seq_along(sets)

    ]


    venn_plot <- VennDiagram::venn.diagram(

        x = sets,

        filename = NULL,

        fill = venn_colors,

        alpha = 0.45,

        cex = 1,

        cat.cex = 0.85,

        cat.fontface = "bold",

        margin = 0.08

    )


    n <- matrix_df$Number_of_Levels


    summary_table <- data.frame(

        Category = c(

            "Common to all 5 levels",

            "Common to exactly 4 levels",

            "Common to exactly 3 levels",

            "Common to exactly 2 levels",

            "Unique to one level",

            "Total unique enriched functions"

        ),

        Number = c(

            sum(n == 5),

            sum(n == 4),

            sum(n == 3),

            sum(n == 2),

            sum(n == 1),

            nrow(matrix_df)

        )

    )


    set_sizes <- sapply(

        expression_levels,

        function(x) {

            sum(
                matrix_df[[x]] == 1
            )

        }

    )


    output_file <- file.path(

        OUTPUT_DIR,

        "Plots",

        paste0(

            database_name,

            "_Venn_and_Summary.png"

        )

    )


    png(

        output_file,

        width = 3200,

        height = 2200,

        res = 250

    )


    grid.newpage()


    # --------------------------------------------------------
    # Title
    # --------------------------------------------------------

    grid.text(

        paste(

            "MIsoTEx - testis",

            database_name

        ),

        x = 0.5,

        y = 0.965,

        gp = gpar(

            fontsize = 23,

            fontface = "bold"

        )

    )


    grid.text(

        paste(

            "Significantly enriched functions | BH adjusted p <",

            ADJ_P_CUTOFF

        ),

        x = 0.5,

        y = 0.935,

        gp = gpar(

            fontsize = 13

        )

    )


    # --------------------------------------------------------
    # Venn
    # --------------------------------------------------------

    pushViewport(

        viewport(

            x = 0.32,

            y = 0.55,

            width = 0.60,

            height = 0.70

        )

    )


    grid.draw(
        venn_plot
    )


    popViewport()


    # --------------------------------------------------------
    # Summary title
    # --------------------------------------------------------

    grid.text(

        "OVERLAP SUMMARY",

        x = 0.82,

        y = 0.82,

        gp = gpar(

            fontsize = 16,

            fontface = "bold"

        )

    )


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary_lines <- paste(

        summary_table$Category,

        ": ",

        summary_table$Number,

        sep = ""

    )


    grid.text(

        paste(

            summary_lines,

            collapse = "\n"

        ),

        x = 0.82,

        y = 0.66,

        gp = gpar(

            fontsize = 11

        )

    )


    # --------------------------------------------------------
    # Set size
    # --------------------------------------------------------

    grid.text(

        "SIGNIFICANT FUNCTIONS PER LEVEL",

        x = 0.82,

        y = 0.43,

        gp = gpar(

            fontsize = 14,

            fontface = "bold"

        )

    )


    set_lines <- paste(

        expression_levels,

        ": ",

        set_sizes,

        sep = ""

    )


    grid.text(

        paste(

            set_lines,

            collapse = "\n"

        ),

        x = 0.82,

        y = 0.30,

        gp = gpar(

            fontsize = 11

        )

    )


    # --------------------------------------------------------
    # Common to all
    # --------------------------------------------------------

    grid.text(

        paste(

            "Functions common to all five levels:",

            sum(n == 5)

        ),

        x = 0.32,

        y = 0.08,

        gp = gpar(

            fontsize = 15,

            fontface = "bold"

        )

    )


    dev.off()


    cat(

        "\nVenn + summary saved:\n",

        output_file,

        "\n"

    )

}


create_venn_summary(

    GO_BP_MATRIX,

    "GO_BP"

)


create_venn_summary(

    GO_MF_MATRIX,

    "GO_MF"

)


create_venn_summary(

    GO_CC_MATRIX,

    "GO_CC"

)


# ============================================================
# 30. GO DOT PLOTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GO DOT PLOTS\n")
cat("============================================================\n")


for (ontology in c(

    "BP",

    "MF",

    "CC"

)) {


    for (level in expression_levels) {


        key <- paste(

            ontology,

            level,

            sep = "_"

        )


        result <- GO_RESULTS[[key]]


        if (is.null(result)) {

            next

        }


        significant <- result %>%

            filter(

                FDR_BH < ADJ_P_CUTOFF

            )


        if (nrow(significant) == 0) {

            next

        }


        genes <- unique(

            transcript_mapping$ENTREZID[

                transcript_mapping$Expression_type ==

                    level &

                !is.na(

                    transcript_mapping$ENTREZID

                )

            ]

        )


        ego_plot <- tryCatch(

            {

                clusterProfiler::enrichGO(

                    gene = genes,

                    universe =
                        background_entrez,

                    OrgDb =
                        org.Hs.eg.db,

                    keyType =
                        "ENTREZID",

                    ont =
                        ontology,

                    pAdjustMethod =
                        "BH",

                    pvalueCutoff =
                        ADJ_P_CUTOFF,

                    qvalueCutoff =
                        1,

                    readable =
                        TRUE

                )

            },

            error = function(e) {

                NULL

            }

        )


        if (is.null(ego_plot)) {

            next

        }


        p <- enrichplot::dotplot(

            ego_plot,

            showCategory =
                TOP_N

        ) +

            ggtitle(

                paste(

                    "MIsoTEx testis - GO",

                    ontology,

                    "-",

                    level

                )

            )


        safe_level <- gsub(

            " ",

            "_",

            level

        )


        ggsave(

            filename = file.path(

                OUTPUT_DIR,

                "Plots",

                paste0(

                    "GO_",

                    ontology,

                    "_",

                    safe_level,

                    "_dotplot.png"

                )

            ),

            plot = p,

            width = 10,

            height = 8,

            dpi = 300

        )

    }

}


# ============================================================
# 31. KEGG ENRICHMENT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("KEGG ENRICHMENT\n")
cat("============================================================\n")


KEGG_RESULTS <- list()


kegg_link <- tryCatch(

    {

        KEGGREST::keggLink(

            "pathway",

            "hsa"

        )

    },

    error = function(e) {

        warning(

            "\nCould not retrieve KEGG information.\n",

            "KEGG analysis will be skipped.\n",

            e$message

        )

        NULL

    }

)


if (!is.null(kegg_link)) {


    kegg_term2gene <- data.frame(

        Pathway = sub(

            "path:",

            "",

            as.character(
                kegg_link
            )

        ),

        Gene = sub(

            "hsa:",

            "",

            names(
                kegg_link
            )

        ),

        stringsAsFactors = FALSE

    )


    pathway_names <- tryCatch(

        {

            KEGGREST::keggList(

                "pathway",

                "hsa"

            )

        },

        error = function(e) {

            NULL

        }

    )


    kegg_term2name <- NULL


    if (!is.null(pathway_names)) {


        kegg_term2name <- data.frame(

            Pathway = sub(

                "path:",

                "",

                names(
                    pathway_names
                )

            ),

            Description = sub(

                " - Homo sapiens \\(human\\)$",

                "",

                as.character(
                    pathway_names
                )

            ),

            stringsAsFactors = FALSE

        )

    }


    write.csv(

        kegg_term2gene,

        file.path(

            OUTPUT_DIR,

            "KEGG",

            "KEGG_TERM2GENE.csv"

        ),

        row.names = FALSE

    )


    if (!is.null(kegg_term2name)) {

        write.csv(

            kegg_term2name,

            file.path(

                OUTPUT_DIR,

                "KEGG",

                "KEGG_TERM2NAME.csv"

            ),

            row.names = FALSE

        )

    }


    # --------------------------------------------------------
    # KEGG per level
    # --------------------------------------------------------

    for (level in expression_levels) {


        level_df <- transcript_mapping %>%

            filter(

                Expression_type == level

            )


        genes <- unique(

            level_df$ENTREZID[

                !is.na(
                    level_df$ENTREZID
                )

            ]

        )


        if (length(genes) == 0) {

            next

        }


        kegg_result <- tryCatch(

            {

                clusterProfiler::enricher(

                    gene = genes,

                    universe =
                        background_entrez,

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

                    "KEGG failed for ",

                    level,

                    ": ",

                    e$message

                )

                NULL

            }

        )


        if (is.null(kegg_result)) {

            next

        }


        result <- as.data.frame(

            kegg_result

        )


        if (nrow(result) == 0) {

            next

        }


        result <- add_statistics(

            result,

            level

        )


        result$Transcript_IDs <- sapply(

            result$geneID,

            get_transcript_contributors,

            level_df = level_df

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


        KEGG_RESULTS[[level]] <-
            result


        safe_level <- gsub(

            " ",

            "_",

            level

        )


        write.csv(

            result,

            file.path(

                OUTPUT_DIR,

                "Individual_KEGG",

                paste0(

                    "KEGG_",

                    safe_level,

                    "_ALL.csv"

                )

            ),

            row.names = FALSE

        )


        write.csv(

            result %>%

                filter(

                    FDR_BH < ADJ_P_CUTOFF

                ),

            file.path(

                OUTPUT_DIR,

                "Individual_KEGG",

                paste0(

                    "KEGG_",

                    safe_level,

                    "_SIGNIFICANT.csv"

                )

            ),

            row.names = FALSE

        )

    }


    if (length(KEGG_RESULTS) > 0) {


        kegg_combined <- bind_rows(

            KEGG_RESULTS

        )


        write.csv(

            kegg_combined,

            file.path(

                OUTPUT_DIR,

                "KEGG",

                "ALL_LEVELS_KEGG.csv"

            ),

            row.names = FALSE

        )


        write.csv(

            kegg_combined %>%

                filter(

                    FDR_BH < ADJ_P_CUTOFF

                ),

            file.path(

                OUTPUT_DIR,

                "KEGG",

                "SIGNIFICANT_ALL_LEVELS_KEGG.csv"

            ),

            row.names = FALSE

        )

    }

}


# ============================================================
# 32. REACTOME ENRICHMENT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("REACTOME ENRICHMENT\n")
cat("============================================================\n")


REACTOME_RESULTS <- list()


for (level in expression_levels) {


    level_df <- transcript_mapping %>%

        filter(

            Expression_type == level

        )


    genes <- unique(

        level_df$ENTREZID[

            !is.na(
                level_df$ENTREZID
            )

        ]

    )


    if (length(genes) == 0) {

        next

    }


    reactome_result <- tryCatch(

        {

            ReactomePA::enrichPathway(

                gene = genes,

                universe =
                    background_entrez,

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

                "Reactome failed for ",

                level,

                ": ",

                e$message

            )

            NULL

        }

    )


    if (is.null(reactome_result)) {

        next

    }


    result <- as.data.frame(

        reactome_result

    )


    if (nrow(result) == 0) {

        next

    }


    result <- add_statistics(

        result,

        level

    )


    result$Transcript_IDs <- sapply(

        result$geneID,

        get_transcript_contributors,

        level_df = level_df

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


    REACTOME_RESULTS[[level]] <-
        result


    safe_level <- gsub(

        " ",

        "_",

        level

    )


    write.csv(

        result,

        file.path(

            OUTPUT_DIR,

            "Individual_Reactome",

            paste0(

                "Reactome_",

                safe_level,

                "_ALL.csv"

            )

        ),

        row.names = FALSE

    )


    write.csv(

        result %>%

            filter(

                FDR_BH < ADJ_P_CUTOFF

            ),

        file.path(

            OUTPUT_DIR,

            "Individual_Reactome",

            paste0(

                "Reactome_",

                safe_level,

                "_SIGNIFICANT.csv"

            )

        ),

        row.names = FALSE

    )

}


if (length(REACTOME_RESULTS) > 0) {


    reactome_combined <- bind_rows(

        REACTOME_RESULTS

    )


    write.csv(

        reactome_combined,

        file.path(

            OUTPUT_DIR,

            "Reactome",

            "ALL_LEVELS_REACTOME.csv"

        ),

        row.names = FALSE

    )


    write.csv(

        reactome_combined %>%

            filter(

                FDR_BH < ADJ_P_CUTOFF

            ),

        file.path(

            OUTPUT_DIR,

            "Reactome",

            "SIGNIFICANT_ALL_LEVELS_REACTOME.csv"

        ),

        row.names = FALSE

    )

}


# ============================================================
# 33. ENRICHMENT SUMMARY TABLE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ENRICHMENT SUMMARY\n")
cat("============================================================\n\n")


go_summary <- data.frame(

    Database = character(0),

    Expression_Type = character(0),

    Total_Terms = integer(0),

    Significant_Terms = integer(0),

    stringsAsFactors = FALSE

)


for (ontology in c(

    "BP",

    "MF",

    "CC"

)) {


    for (level in expression_levels) {


        result <- GO_RESULTS[[

            paste(

                ontology,

                level,

                sep = "_"

            )

        ]]


        if (is.null(result)) {

            total <- 0

            significant <- 0

        } else {

            total <- nrow(
                result
            )

            significant <- sum(

                result$FDR_BH <

                    ADJ_P_CUTOFF,

                na.rm = TRUE

            )

        }


        go_summary <- rbind(

            go_summary,

            data.frame(

                Database =

                    paste0(

                        "GO-",

                        ontology

                    ),

                Expression_Type =

                    level,

                Total_Terms =

                    total,

                Significant_Terms =

                    significant,

                stringsAsFactors =

                    FALSE

            )

        )

    }

}


print(
    go_summary
)


write.csv(

    go_summary,

    file.path(

        OUTPUT_DIR,

        "Comparison",

        "GO_Enrichment_Summary.csv"

    ),

    row.names = FALSE

)


# ============================================================
# 34. FINAL OUTPUT SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("MIsoTEx FUNCTIONAL ENRICHMENT ANALYSIS COMPLETE\n")
cat("============================================================\n\n")


cat(

    "Input:\n",

    INPUT_FILE,

    "\n\n"

)


cat(

    "Total transcripts analyzed: ",

    nrow(df),

    "\n"

)


cat(

    "No >50 transcript/gene filter applied.\n"

)


cat(

    "Enrichment unit: unique Entrez genes.\n"

)


cat(

    "Background: complete human Entrez gene set.\n"

)


cat(

    "Multiple-testing correction: Benjamini-Hochberg.\n"

)


cat(

    "Significance cutoff: FDR_BH < ",

    ADJ_P_CUTOFF,

    "\n\n"

)


cat(

    "Output directory:\n",

    normalizePath(
        OUTPUT_DIR
    ),

    "\n\n"

)


cat("IMPORTANT OUTPUTS:\n\n")


cat(

    "1. GO-BP UpSet:\n",

    file.path(

        OUTPUT_DIR,

        "Plots",

        "GO_BP_UpSet_Plot.png"

    ),

    "\n\n"

)


cat(

    "2. GO-BP Venn:\n",

    file.path(

        OUTPUT_DIR,

        "Plots",

        "GO_BP_Venn.png"

    ),

    "\n\n"

)


cat(

    "3. GO-BP Venn + Summary:\n",

    file.path(

        OUTPUT_DIR,

        "Plots",

        "GO_BP_Venn_and_Summary.png"

    ),

    "\n\n"

)


cat(

    "4. GO-BP function matrix:\n",

    file.path(

        OUTPUT_DIR,

        "Comparison",

        "GO_BP_Function_Presence_Absence_Matrix.csv"

    ),

    "\n\n"

)


cat(

    "5. GO-BP intersection table:\n",

    file.path(

        OUTPUT_DIR,

        "Comparison",

        "GO_BP_UpSet_Intersection_Summary.csv"

    ),

    "\n\n"

)


cat(

    "6. GO-BP overlap summary:\n",

    file.path(

        OUTPUT_DIR,

        "Comparison",

        "GO_BP_Overlap_Summary.csv"

    ),

    "\n\n"

)


cat(

    "7. GO-BP significant enrichment:\n",

    file.path(

        OUTPUT_DIR,

        "GO_BP",

        "SIGNIFICANT_ALL_LEVELS.csv"

    ),

    "\n\n"

)


cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")
