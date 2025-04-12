# --------------------------------------------------------------------------
# ZADANIE: Analiza Różnicowej Ekspresji Genów (DEG) i Analiza Funkcjonalna
#          dla danych RNA-Seq ludzkiego mitochondrium
# --------------------------------------------------------------------------
# Cel: Przeprowadzenie pełnej analizy DEG przy użyciu pakietu DESeq2
#      oraz analizy funkcjonalnej z gprofiler2 na danych zliczeniowych
#      uzyskanych na poprzednich zajęciach (wynik featureCounts).
#      Skrypt bazuje na przykładzie analizy danych 'pasilla'.
# --------------------------------------------------------------------------

# 1. ŁADOWANIE WYMAGANYCH PAKIETÓW
# --------------------------------------------------------------------------
# Upewnij się, że wszystkie pakiety są zainstalowane.
# (Kod instalacyjny można skopiować z poprzedniego skryptu lub założyć, że są zainstalowane)

setwd("") # Ustaw własną ścieżkę do katalogu roboczego

library("DESeq2")       # Analiza różnicowej ekspresji
library("ggplot2")      # Zaawansowane wizualizacje
library("pheatmap")     # Tworzenie map ciepła (heatmaps)
library("RColorBrewer") # Palety kolorów do wizualizacji
library("dplyr")        # Manipulacja danymi
library("gprofiler2")   # Analiza funkcjonalna GO i KEGG
library("Rsubread")       # 
# --------------------------------------------------------------------------
# 2. PRZYGOTOWANIE DANYCH WEJŚCIOWYCH
# --------------------------------------------------------------------------

buildindex(basename="HumanMT", 
           reference="Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz")

align.stat_113 <- align(index="HumanMT",readfile1="reads/SRR2816113.fastq.gz",
                        output_file="HsMT_113.BAM",phredOffset=64)

align.stat_114 <- align(index="HumanMT",readfile1="reads/SRR2816114.fastq.gz",
                        output_file="HsMT_114.BAM",phredOffset=64)

align.stat_115 <- align(index="HumanMT",readfile1="reads/SRR2816115.fastq.gz",
                        output_file="HsMT_115.BAM",phredOffset=64)

align.stat_113



#Zliczanie countów ----
ann <- "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gtf.gz"
fc <- featureCounts(files = c("HsMT_113.BAM",
                              "HsMT_114.BAM",
                              "HsMT_115.BAM"),
                    annot.ext = "Homo_sapiens.GRCh38.113.chromosome.MT.gff3",
                    isGTFAnnotationFile = TRUE,
                    GTF.featureType = "gene",
                    #GTF.attrType = "Parent",    
                    useMetaFeatures = TRUE,
                    nthreads = 4)

colnames(fc$counts) <- gsub(x = colnames(fc$counts), pattern = "\\.BAM",  replacement = "")
fc$counts           
head(fc$counts)     
fc$annotation
head(fc)


# Przypisanie macierzy zliczeń do zmiennej 'cts' dla spójności z przykładem pasilla
cts <- fc$counts

# Wyświetlenie wymiarów i początku macierzy zliczeń
print(dim(cts))
print(head(cts))

# --- Tworzenie tabeli metadanych (colData) ---

# !!! WAŻNE !!!
# Poniższa tabela 'coldata' jest PRZYKŁADOWA.
# Należy ją dostosować do RZECZYWISTEGO planu eksperymentu dla próbek
# SRR2816113, SRR2816114, SRR2816115.

# Założenie przykładowe:
# HsMT_113: Warunek A (np. kontrola)
# HsMT_114: Warunek B (np. traktowanie 1)
# HsMT_115: Warunek B (np. traktowanie 1 - replika)

# Pobranie nazw próbek z macierzy zliczeń
sample_names <- colnames(cts)
print(paste("Nazwy próbek:", paste(sample_names, collapse=", ")))

# Stworzenie przykładowej ramki danych colData
# Upewnij się, że nazwy wierszy w coldata DOKŁADNIE odpowiadają nazwom kolumn w cts!
coldata <- data.frame(
  condition = factor(c("treated", "untreated", "untreated")) # PRZYKŁADOWE warunki - dostosuj!
  # Można dodać inne kolumny, np. typ biblioteki, batch itp.
  # type = factor(c("Type1", "Type1", "Type1"))
)
rownames(coldata) <- sample_names

# Sprawdzenie zgodności nazw wierszy coldata i kolumn cts
if (!all(rownames(coldata) == colnames(cts))) {
  stop("BŁĄD: Nazwy wierszy w 'coldata' nie pasują do nazw kolumn w 'cts'!")
}

# Wyświetlenie tabeli metadanych
print(coldata)

# --------------------------------------------------------------------------
# 3. TWORZENIE OBIEKTU DESeqDataSet
# --------------------------------------------------------------------------
# Tworzymy obiekt DESeqDataSet, łącząc macierz zliczeń (cts),
# tabelę metadanych (coldata) i formułę projektu (design).
# Formuła `~ condition` wskazuje, że chcemy porównać grupy zdefiniowane
# w kolumnie 'condition' tabeli coldata.

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition) # Używamy kolumny 'condition' zdefiniowanej w coldata

# Wyświetlenie podsumowania obiektu DESeqDataSet
print(dds)

# --------------------------------------------------------------------------
# 4. USTAWIANIE POZIOMU REFERENCYJNEGO (ZALECANE)
# --------------------------------------------------------------------------
# Ustawiamy poziom referencyjny dla porównań. Jeśli porównujemy CondA vs CondB,
# ustawienie CondA jako referencji sprawi, że log2FoldChange będzie oznaczał
# zmianę w CondB względem CondA.

# Ustawienie poziomu referencyjnego (dostosuj do swoich warunków!)
# Np. jeśli chcesz porównać CondB vs CondA:
dds$condition <- relevel(dds$condition, ref = "untreated")

# Sprawdzenie poziomów faktora
print(levels(dds$condition))

# --------------------------------------------------------------------------
# 5. URUCHOMIENIE ANALIZY DESeq2
# --------------------------------------------------------------------------
# Uruchomienie głównej funkcji DESeq(), która wykonuje normalizację,
# estymację dyspersji i testowanie statystyczne.

dds <- DESeq(dds)

# --------------------------------------------------------------------------
# 6. POBRANIE I FILTROWANIE WYNIKÓW
# --------------------------------------------------------------------------
# Pobranie wyników dla domyślnego porównania (zdefiniowanego przez relevel)
res <- results(dds)

# Wyświetlenie podsumowania wyników
summary(res)

# Konwersja na data.frame i dodanie ID genów (Ensembl IDs z nazw wierszy)
res_tab <- as.data.frame(res)
res_tab$gene_id <- rownames(res_tab)

# Usunięcie wierszy z NA w padj
res_sig <- res_tab[!is.na(res_tab$padj),]
res_sig
# Definicja progów istotności
log2FC_threshold <- 0.1 # Próg dla log2 krotności zmiany (np. 1.0 dla 2-krotnej zmiany)
padj_threshold <- 0.5   # Próg dla skorygowanej wartości p (FDR)

# Identyfikacja genów różnicowo wyeksponowanych (DEG)
DEG <- res_sig[abs(res_sig$log2FoldChange) > log2FC_threshold & res_sig$padj < padj_threshold, ]

# Wyświetlenie liczby znalezionych DEG i ich początku
cat("\nLiczba zidentyfikowanych genów DEG:", nrow(DEG), "\n")
print(head(DEG))

# --------------------------------------------------------------------------
# 7. WIZUALIZACJA WYNIKÓW
# --------------------------------------------------------------------------
# Wizualizacje pomagają zrozumieć wyniki i ocenić jakość analizy.

# --- 7.1. MA Plot ---
plotMA(res, ylim=c(-5,5)) # Dostosuj ylim w razie potrzeby
title("MA Plot")

# --- 7.2. Dispersion Plot ---
plotDispEsts(dds)
title("Dispersion Estimates Plot")

# --- 7.3. Transformacja danych (np. VST) ---
# Potrzebna do PCA i heatmap
vsd <- vst(dds, blind=FALSE)

# --- 7.4. PCA Plot ---
# Używamy 'condition' do kolorowania punktów (lub innej zmiennej z coldata)
pca_plot <- plotPCA(vsd, intgroup="condition")
pca_plot <- pca_plot + ggtitle("PCA Plot (VST Data)") + theme_bw()
print(pca_plot)

# --- 7.5. Sample Distance Heatmap ---
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, rownames(vsd@colData), sep="-") # Lepsze etykiety
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Sample Distance Heatmap (VST Data)")

# --- 7.6. Volcano Plot ---
# Przygotowanie danych (dodanie kolumny 'significant')
res_sig <- res_sig %>%
  mutate(significant = ifelse(padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold,
                              ifelse(log2FoldChange > log2FC_threshold, "Up-regulated", "Down-regulated"),
                              "Not Significant"))

# Tworzenie wykresu
volcano_plot <- ggplot(res_sig, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_bw() + theme(legend.position = "bottom")
print(volcano_plot)

# --- 7.7. Counts Plot for Single Gene ---
# Wykres dla genu o najwyższej zmianie (jeśli istnieją DEG)
if (nrow(DEG) > 0) {
  top_gene <- DEG[order(-abs(DEG$log2FoldChange)), "gene_id"][1] # Gen o największej bezwzględnej zmianie
  gene_counts <- plotCounts(dds, gene = top_gene, intgroup = "condition", returnData = TRUE)
  counts_plot <- ggplot(gene_counts, aes(x = condition, y = count, color = condition)) + # Kolor wg condition
    geom_point(position = position_jitter(width = 0.1), size = 3) +
    scale_y_log10() +
    labs(title = paste("Normalized Counts for Gene:", top_gene), y = "Normalized Counts (log10 scale)") +
    theme_bw() +
    theme(legend.position = "none") # Ukrycie legendy dla kolorów, jeśli zbędna
  print(counts_plot)
} else {
  print("Nie znaleziono genów DEG, pominięto wykres zliczeń.")
}

# --------------------------------------------------------------------------
# 8. ANALIZA FUNKCJONALNA (GO i KEGG) z gProfiler2
# --------------------------------------------------------------------------
# Przeprowadzamy analizę funkcjonalną dla listy zidentyfikowanych genów DEG.

if (nrow(DEG) > 0) {
  # Wyodrębnienie identyfikatorów genów DEG (zakładamy, że są to Ensembl IDs)
  deg_gene_ids <- DEG$gene_id
  
  # Uruchomienie analizy funkcjonalnej dla Homo sapiens
  gost_results <- gost(query = deg_gene_ids,
                       organism = "hsapiens", # !!! Ważne: organizm to człowiek !!!
                       sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"), # Szerszy zakres baz
                       evcodes = TRUE,
                       user_threshold = 0.05,
                       correction_method = "g_SCS")
  
  # Sprawdzenie i wyświetlenie wyników
  if (!is.null(gost_results) && nrow(gost_results$result) > 0) {
    functional_analysis_results <- gost_results$result
    print(head(functional_analysis_results[order(functional_analysis_results$p_value),
                                           c("source", "term_id", "term_name", "p_value", "intersection_size")]))
    
    # Wizualizacja wyników (opcjonalnie)
    gost_plot_static <- gostplot(gost_results, capped = TRUE, interactive = FALSE)
    print(gost_plot_static)
    
  } else {
    print("Nie znaleziono istotnych wyników analizy funkcjonalnej.")
  }
  
} else {
  print("Pominięto analizę funkcjonalną, ponieważ nie znaleziono genów DEG.")
}

# --------------------------------------------------------------------------
# KONIEC ZADANIA
# --------------------------------------------------------------------------
# Zapisz wyniki, wykresy lub cały skrypt.
# Np. zapis tabeli DEG do pliku CSV:
# write.csv(DEG, file = "DEG_results.csv", row.names = FALSE)
# --------------------------------------------------------------------------

