# --------------------------------------------------------------------------
# Analiza różnicowej ekspresji genów (DEG) przy użyciu DESeq2 na danych pasilla
# --------------------------------------------------------------------------

# 1. INSTALACJA I ŁADOWANIE WYMAGANYCH PAKIETÓW
# --------------------------------------------------------------------------

# Sprawdzenie i instalacja BiocManager, jeśli nie jest zainstalowany
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalacja pakietu 'pasilla' zawierającego dane eksperymentalne, jeśli nie jest zainstalowany
# Pakiet 'pasilla' dostarcza przykładowy zestaw danych RNA-Seq z eksperymentu na muszkach owocowych Drosophila melanogaster.
if (!requireNamespace("pasilla", quietly = TRUE))
  BiocManager::install("pasilla")

# Instalacja pakietu 'DESeq2' do analizy różnicowej ekspresji genów, jeśli nie jest zainstalowany
# DESeq2 jest standardowym narzędziem do normalizacji i analizy danych zliczeń RNA-Seq.
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

# Instalacja pakietów do wizualizacji, jeśli nie są zainstalowane
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  install.packages("RColorBrewer")
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) # Dodano dplyr (użyty w volcano plot)
  install.packages("dplyr")
if (!requireNamespace("gprofiler2", quietly = TRUE)) # Dodano gprofiler2
  install.packages("gprofiler2")


# Ładowanie zainstalowanych bibliotek do bieżącej sesji R
library(dplyr)
library("pasilla")      # Dane eksperymentalne
library("DESeq2")       # Analiza różnicowej ekspresji
library("RColorBrewer") # Palety kolorów do wizualizacji
library("pheatmap")     # Tworzenie map ciepła (heatmaps)
library("ggplot2")      # Zaawansowane wizualizacje (np. volcano plot)
library("dplyr")        # Manipulacja danymi (użyty w volcano plot)
library("gprofiler2")   # Analiza wzbogacenia funkcjonalnego GO i KEGG

# --------------------------------------------------------------------------
# 2. WCZYTYWANIE I PRZYGOTOWANIE DANYCH
# --------------------------------------------------------------------------

# Wczytywanie ścieżek do plików z danymi zliczeń i adnotacjami próbek z pakietu 'pasilla'
# system.file() lokalizuje pliki w zainstalowanych pakietach R.
# mustWork=TRUE zapewnia, że funkcja zwróci błąd, jeśli plik nie zostanie znaleziony.

# Ścieżka do pliku z surowymi danymi zliczeń genów (raw read counts)
pasCtsFile <- system.file("extdata",
                          "pasilla_gene_counts.tsv",
                          package="pasilla", mustWork=TRUE)

# Ścieżka do pliku z adnotacjami próbek (metadane eksperymentu)
pasAnnoFile <- system.file("extdata",
                           "pasilla_sample_annotation.csv",
                           package="pasilla", mustWork=TRUE)

# Wczytywanie danych zliczeń genów z pliku TSV (tab-separated values)
# sep="\t" określa tabulator jako separator kolumn.
# row.names="gene_id" ustawia pierwszą kolumnę (identyfikatory genów) jako nazwy wierszy.
# as.matrix() konwertuje wczytaną ramkę danych (data frame) na macierz, co jest wymagane przez DESeq2.
cts <- as.matrix(read.csv(pasCtsFile, sep="\t", row.names="gene_id"))

# Wczytywanie adnotacji próbek (metadanych) z pliku CSV
# row.names=1 ustawia pierwszą kolumnę (identyfikatory próbek) jako nazwy wierszy.
coldata_raw <- read.csv(pasAnnoFile, row.names=1)

# Wybór interesujących nas kolumn z metadanych: 'condition' (warunek eksperymentalny) i 'type' (typ biblioteki)
coldata <- coldata_raw[, c("condition", "type")]

# Konwersja kolumn 'condition' i 'type' na typ 'factor' (kategoryczny)
# Faktory są używane w R do reprezentowania zmiennych kategorycznych i są ważne dla modelowania statystycznego w DESeq2.
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# Wyświetlenie pierwszych 10 wierszy macierzy zliczeń (dla wglądu w dane)
head(cts, 10)

# Wyświetlenie przygotowanej tabeli z metadanymi (dla wglądu)
print(coldata)

# --------------------------------------------------------------------------
# 3. SPRAWDZENIE I DOPASOWANIE NAZW PRÓBEK
# --------------------------------------------------------------------------

# W oryginalnych danych nazwy wierszy w 'coldata' (metadane) mogą mieć prefiks 'fb',
# podczas gdy nazwy kolumn w 'cts' (zliczenia) go nie mają. Musimy je ujednolicić.

# Usunięcie prefiksu "fb" z nazw wierszy w tabeli metadanych 'coldata'
rownames(coldata) <- sub("fb", "", rownames(coldata))

# Sprawdzenie, czy wszystkie nazwy próbek z metadanych ('coldata') znajdują się wśród nazw kolumn macierzy zliczeń ('cts')
# Powinno zwrócić TRUE.
all(rownames(coldata) %in% colnames(cts))

# Sprawdzenie, czy nazwy próbek w 'coldata' i 'cts' są DOKŁADNIE w tej samej kolejności.
# Prawdopodobnie zwróci FALSE, co oznacza, że musimy uporządkować kolumny macierzy zliczeń.
all(rownames(coldata) == colnames(cts))

# Uporządkowanie kolumn macierzy zliczeń ('cts') tak, aby odpowiadały kolejności wierszy w tabeli metadanych ('coldata')
# Jest to kluczowe dla poprawnego działania DESeq2, który oczekuje zgodności między zliczeniami a metadanymi.
cts <- cts[, rownames(coldata)]

# Ponowne sprawdzenie, czy nazwy próbek są teraz w tej samej kolejności.
# Powinno zwrócić TRUE.
all(rownames(coldata) == colnames(cts))

# --------------------------------------------------------------------------
# 4. TWORZENIE OBIEKTU DESeqDataSet
# --------------------------------------------------------------------------

# Tworzenie obiektu DESeqDataSet, który jest podstawową strukturą danych używaną przez pakiet DESeq2.
# Łączy on macierz zliczeń ('countData'), tabelę metadanych ('colData') oraz formułę projektową ('design').
# Formuła projektowa (`~ condition`) określa, które zmienne z metadanych będą używane do modelowania ekspresji genów.
# W tym przypadku modelujemy wpływ zmiennej 'condition' (np. treated vs untreated).
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Wyświetlenie podsumowania obiektu DESeqDataSet
print(dds)

# --------------------------------------------------------------------------
# 5. USTAWIANIE POZIOMU REFERENCYJNEGO
# --------------------------------------------------------------------------

# Upewnienie się, że poziom 'untreated' jest poziomem referencyjnym dla faktora 'condition'.
# DESeq2 domyślnie wybiera poziom referencyjny alfabetycznie. Jawne ustawienie referencji
# zapewnia, że wyniki log2FoldChange będą interpretowane jako zmiana WZGLĘDEM grupy 'untreated'.

# Można to zrobić na dwa sposoby:

# Sposób 1: Ustawienie poziomów faktora w określonej kolejności
# dds$condition <- factor(dds$condition, levels = c("untreated", "treated"))

# Sposób 2: Użycie funkcji relevel() do jawnego ustawienia poziomu referencyjnego
dds$condition <- relevel(dds$condition, ref = "untreated")

# Sprawdzenie, czy poziom referencyjny został ustawiony poprawnie (pierwszy poziom powinien być "untreated")
levels(dds$condition)

# --------------------------------------------------------------------------
# 6. URUCHOMIENIE ANALIZY DESeq2
# --------------------------------------------------------------------------

# Uruchomienie głównej funkcji DESeq(), która wykonuje następujące kroki:
# 1. Estymacja czynników normalizacyjnych (size factors) w celu uwzględnienia różnic w głębokości sekwencjonowania.
# 2. Estymacja dyspersji genów (gene-wise dispersion estimates).
# 3. Dopasowanie uogólnionego modelu liniowego (GLM) i testowanie istotności statystycznej (Wald test).
dds <- DESeq(dds)

# Pobranie wyników analizy różnicowej ekspresji dla porównania zdefiniowanego przez projekt (domyślnie porównuje poziomy faktora 'condition').
# Wyniki zawierają m.in. log2 krotności zmiany (log2FoldChange), wartość p (pvalue) i skorygowaną wartość p (padj).
res <- results(dds)

# Wyświetlenie pierwszych kilku wierszy tabeli wyników
head(res)

# Wyświetlenie całkowitej liczby genów w wynikach (przed filtrowaniem)
nrow(res)

# --------------------------------------------------------------------------
# 7. FILTROWANIE I EKSPLORACJA WYNIKÓW
# --------------------------------------------------------------------------

# Konwersja obiektu wyników DESeqResults na standardową ramkę danych (data frame) dla łatwiejszej manipulacji
res_tab <- as.data.frame(res)

# Dodanie kolumny z nazwami genów (identyfikatorami) do ramki danych wyników
res_tab$gene_id <- rownames(res_tab)

# Filtrowanie wyników:
# 1. Usunięcie wierszy (genów), dla których nie udało się obliczyć log2FoldChange (np. z powodu zerowych zliczeń we wszystkich próbkach jednej grupy).
res_sig <- res_tab[!is.na(res_tab$log2FoldChange),]
# 2. Usunięcie wierszy, dla których nie udało się obliczyć skorygowanej wartości p (padj) (często z powodu niskiej liczby zliczeń lub dużej dyspersji).
res_sig <- res_sig[!is.na(res_sig$padj),]

# Identyfikacja genów różnicowo wyeksponowanych (Differentially Expressed Genes - DEG)
# Kryteria:
# - Bezwzględna wartość log2 krotności zmiany (|log2FoldChange|) > 1 (co odpowiada co najmniej 2-krotnej zmianie ekspresji)
# - Skorygowana wartość p (padj) < 0.05 (kontrola fałszywych odkryć - FDR)
DEG <- res_sig[abs(res_sig$log2FoldChange) > 1 & res_sig$padj < 0.05, ]

# Wyświetlenie liczby zidentyfikowanych genów DEG
nrow(DEG)

# Wyświetlenie pierwszych kilku wierszy tabeli z genami DEG
head(DEG)

# Przykład: Wybór genów z log2FoldChange w określonym zakresie (np. blisko 1)
selected <- res_sig[res_sig$log2FoldChange > 0.9 & res_sig$log2FoldChange < 1.1,]

# Wyświetlenie surowych zliczeń dla genów wybranych w poprzednim kroku
# Sprawdzenie, które wiersze macierzy 'cts' mają nazwy odpowiadające genom w 'selected'
cts[rownames(cts) %in% rownames(selected),]

# Przykład: Wybór konkretnych genów na podstawie ich identyfikatorów FlyBase (FBgn)
sel <- c("FBgn0000078", "FBgn0000079", "FBgn0000233", "FBgn0000359")

# Wyświetlenie surowych zliczeń dla wybranych genów
cts[rownames(cts) %in% sel,]

# Wyświetlenie wyników DESeq2 dla tych konkretnych genów
res_sig[rownames(res_sig) %in% sel,]

# Przykład obliczenia logarytmu o podstawie 2
log2(0.5) # Oczekiwany wynik: -1

# --------------------------------------------------------------------------
# 8. WIZUALIZACJA WYNIKÓW
# --------------------------------------------------------------------------

# --- 8.1. Wykres MA (MA plot) ---
# Wizualizuje zależność log2 krotności zmiany (M) od średniej znormalizowanej ekspresji (A).
# Pokazuje rozkład zmian ekspresji genów i wyróżnia geny istotnie różnicowo wyeksponowane (domyślnie na czerwono).
# ylim=c(-2,2) ogranicza oś Y dla lepszej widoczności większości genów.
plotMA(res, ylim=c(-2,2))
title("MA Plot") # Dodanie tytułu

# --- 8.2. Wykres rozrzutu dyspersji (Dispersion Plot) ---
# Pokazuje oszacowania dyspersji dla każdego genu (czarne kropki), dopasowaną krzywą zależności dyspersji od średniej (czerwona linia)
# oraz ostateczne oszacowania dyspersji użyte w modelu (niebieskie kropki, jeśli były przesunięte w kierunku krzywej).
# Jest to ważny wykres diagnostyczny pokazujący, jak model radzi sobie ze zmiennością danych.
plotDispEsts(dds)
title("Dispersion Estimates Plot") # Dodanie tytułu


# --- 8.3. Transformacja danych ---
# Surowe dane zliczeń mają rozkład, który nie jest idealny do niektórych wizualizacji (np. PCA, heatmapy).
# Transformacje stabilizujące wariancję (VST) lub logarytmiczna (rlog) pomagają znormalizować dane.

# Variance Stabilizing Transformation (VST) - szybsza dla dużych zbiorów danych
vsd <- vst(dds, blind=FALSE) # blind=FALSE oznacza, że transformacja uwzględnia projekt eksperymentu

# Regularized Logarithm Transformation (rlog) - może być lepsza dla małych zbiorów danych
# rld <- rlog(dds, blind=FALSE) # Odkomentuj, jeśli chcesz użyć rlog

# --- 8.4. Analiza głównych składowych (PCA plot ang. Principal Component Analysis) ---
# PCA redukuje wymiarowość danych, pozwalając zwizualizować podobieństwo między próbkami w 2D.
# Oczekujemy, że próbki z tej samej grupy eksperymentalnej (np. 'condition') będą grupować się razem.
# intgroup określa, które zmienne z metadanych użyć do kolorowania/kształtowania punktów na wykresie.
pca_plot <- plotPCA(vsd, intgroup=c("condition", "type"))
pca_plot <- pca_plot + ggtitle("PCA Plot (VST Data)") + theme_bw() # Dodanie tytułu i poprawa wyglądu
print(pca_plot) # Wyświetlenie zmodyfikowanego wykresu
# Jeśli użyłeś rlog: plotPCA(rld, intgroup=c("condition", "type"))

# --- 8.5. Mapa ciepła odległości między próbkami (Sample Distance Heatmap) ---
# Wizualizuje macierz odległości euklidesowych między próbkami na podstawie przetransformowanych danych ekspresji.
# Podobne próbki będą miały mniejszą odległość (ciemniejszy kolor na heatmapie).

# Obliczenie macierzy odległości euklidesowych między próbkami (transponujemy macierz assay(vsd), aby próbki były w wierszach)
sampleDists <- dist(t(assay(vsd)))

# Konwersja obiektu 'dist' na standardową macierz
sampleDistMatrix <- as.matrix(sampleDists)

# Ustawienie nazw wierszy macierzy odległości, aby zawierały informacje o warunku i typie próbki
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
# Usunięcie nazw kolumn dla przejrzystości (są takie same jak wierszy)
colnames(sampleDistMatrix) <- NULL

# Definicja palety kolorów (od białego do niebieskiego)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Wygenerowanie heatmapy
# clustering_distance_rows/cols - używamy obliczonych odległości do grupowania (klastrowania) wierszy i kolumn
# col - definiuje użytą paletę kolorów
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Sample Distance Heatmap (VST Data)") # Dodanie tytułu

# --- 8.6. heatmapa ekspresji wybranych genów (Gene Expression Heatmap) ---
# Wizualizuje poziomy ekspresji dla wybranej grupy genów we wszystkich próbkach.

# Normalizacja danych zliczeń (alternatywna do VST/rlog, często używana do heatmap genów)
ntd <- normTransform(dds) # Log2(n + 1) transformacja znormalizowanych zliczeń

# Wybór np. 20 genów o najwyższej średniej znormalizowanej ekspresji
# counts(dds, normalized=TRUE) pobiera znormalizowane zliczenia
# rowMeans oblicza średnią ekspresję dla każdego genu
# order(..., decreasing=TRUE) sortuje geny od najwyższej do najniższej średniej ekspresji
# [1:20] wybiera indeksy pierwszych 20 genów
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

# Przygotowanie ramki danych z adnotacjami dla kolumn heatmapy (próbek)
df <- as.data.frame(colData(dds)[,c("condition","type")])

# Wygenerowanie heatmapy dla wybranych genów
# assay(ntd)[select,] wybiera wiersze (geny) z macierzy znormalizowanych danych
# cluster_rows=FALSE - nie grupuj genów (są już posortowane wg ekspresji)
# show_rownames=FALSE - nie pokazuj nazw genów (może być ich dużo)
# cluster_cols=FALSE - nie grupuj próbek (zachowaj oryginalną kolejność lub uporządkuj wg metadanych)
# annotation_col=df - dodaj adnotacje (paski kolorów) dla próbek na podstawie metadanych
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df,
         main="Top 20 Expressed Genes Heatmap (Normalized Counts)") # Dodanie tytułu

# --- 8.7. Volcano Plot ---
# Używa wyników z ramki danych 'res_sig' (po usunięciu NA)

# Definicja progów
log2FC_threshold <- 1.0
padj_threshold <- 0.05

# Dodanie kolumny wskazującej, czy gen ulega istotnej ekspresji
# Używamy funkcji mutate z pakietu dplyr
res_sig <- res_sig %>%
  mutate(significant = ifelse(padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold,
                              ifelse(log2FoldChange > log2FC_threshold, "Up-regulated", "Down-regulated"),
                              "Not Significant"))

# Tworzenie wykresu Volcano Plot za pomocą ggplot2
volcano_plot <- ggplot(res_sig, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) + # Punkty reprezentujące geny
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")) + # Kolory dla grup
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "grey40") + # Linie pionowe dla progu log2FC
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey40") + # Linia pozioma dla progu padj
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Significance") + # Etykiety i legenda
  theme_bw() + # Ustawienie stylu wykresu
  theme(legend.position = "bottom") # Pozycja legendy

print(volcano_plot) # Wyświetlenie wykresu

# --- 8.8. Wykres zliczeń dla pojedynczego genu ---
# Wybierzmy gen o największej zmianie ekspresji (najbardziej up-regulowany)

# Sprawdzenie, czy istnieją jakiekolwiek geny DEG, zanim wybierzemy najlepszy
if (nrow(DEG) > 0) {
  # Sortowanie DEG malejąco wg log2FoldChange i wybranie pierwszego genu
  top_gene <- DEG[order(-DEG$log2FoldChange), "gene_id"][1] # Używamy kolumny gene_id dodanej wcześniej
  
  # Pobranie znormalizowanych danych zliczeń dla tego genu
  gene_counts <- plotCounts(dds, gene = top_gene, intgroup = c("condition", "type"), returnData = TRUE)
  
  # Tworzenie wykresu zliczeń za pomocą ggplot2
  counts_plot <- ggplot(gene_counts, aes(x = condition, y = count, color = type)) +
    geom_point(position = position_jitter(width = 0.1), size = 3) + # Punkty dla każdej próbki, z lekkim rozrzutem
    scale_y_log10() + # Skala logarytmiczna dla osi Y (często lepsza dla zliczeń) - usunięto ręczne breaki dla ogólności
    labs(title = paste("Normalized Counts for Gene:", top_gene),
         x = "Condition",
         y = "Normalized Counts (log10 scale)",
         color = "Library Type") +
    theme_bw()
  
  print(counts_plot) # Wyświetlenie wykresu
} else {
  print("No differentially expressed genes found based on the criteria.")
}


# --------------------------------------------------------------------------
# 9. Adnotacja funkcjonalna (GO i KEGG)
# --------------------------------------------------------------------------

# Sprawdzenie, czy zidentyfikowano jakiekolwiek geny DEG
if (nrow(DEG) > 0) {
  # Przygotowanie listy identyfikatorów genów DEG dla gProfiler
  deg_gene_ids <- DEG$gene_id
  
  # Uruchomienie adnotacji za pomocą funkcji gost() z gprofiler2
  # query: wektor identyfikatorów genów DEG
  # organism: organizm, dla którego przeprowadzana jest analiza (Drosophila melanogaster)
  # sources: bazy danych do przeszukania ('GO:BP' - Biological Process, 'GO:MF' - Molecular Function, 'GO:CC' - Cellular Component, 'KEGG' - ścieżki metaboliczne i sygnałowe)
  # evcodes: czy dołączyć kody dowodowe dla adnotacji GO
  # user_threshold: próg istotności (skorygowanej wartości p)
  # correction_method: metoda korekcji testów wielokrotnych (domyślna 'g_SCS' jest zalecana)
  gost_results <- gost(query = deg_gene_ids,
                       organism = "dmelanogaster",
                       sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
                       evcodes = TRUE,
                       user_threshold = 0.05,
                       correction_method = "g_SCS")
  
  # Sprawdzenie, czy uzyskano wyniki
  if (!is.null(gost_results) && nrow(gost_results$result) > 0) {
    # Wyodrębnienie tabeli wyników
    enrichment_results <- gost_results$result
    
    # Wyświetlenie pierwszych kilku wierszy istotnych wyników adnotacji
    # Sortowanie wg skorygowanej wartości p (p_value)
    head(enrichment_results[order(enrichment_results$p_value),
                            c("query", "source", "term_id", "term_name", "p_value", "intersection_size")])
    
    # Opcjonalnie: Wizualizacja wyników za pomocą gostplot
    # Tworzy interaktywny wykres słupkowy pokazujący zadnotowane funkcje
    # (może wymagać zainstalowania pakietu 'plotly')
    # gost_plot <- gostplot(gost_results, capped = TRUE, interactive = TRUE) # Ustaw interactive = FALSE dla statycznego wykresu ggplot
    # print(gost_plot)
    
    # Generowanie statycznego wykresu ggplot (bardziej niezawodne w różnych środowiskach)
    gost_plot_static <- gostplot(gost_results, capped = TRUE, interactive = TRUE)
    print(gost_plot_static)
    
    
  } else {
    print("No significant enrichment results found for the given DEG list and thresholds.")
  }
  
} else {
  print("Skipping enrichment analysis because no differentially expressed genes were found.")
}


# --------------------------------------------------------------------------
# KONIEC SKRYPTU
# --------------------------------------------------------------------------
