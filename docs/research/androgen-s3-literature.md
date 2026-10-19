# Androgen program, S3 identity and IRI: literature review (scout, 2026-10-18)

A review for notebook 111. Tags: [full text] = methods, results and figures read; [abstract only] = no method or
result claims are taken from it. The scout did not look at any IR effect in our data, only at gene availability. The
two external datasets it found are now in `data/external/xiong2023_ar_ko/` and
`data/external/vinas2020_gse144622/`.

## Verdict

The hypothesis is half supported, and its "rather than" is a false dichotomy.
1. **AR-dependent genes are S2/S3-weighted and S3-maximal.** Sex-biased genes in general are not clearly S3-heaviest:
   two bulk sources put S2 at or above S3.
2. **The male PT program collapses after IRI.** This is supported, including in published but unremarked data
   (Kirita's AR regulon).
3. **It is not "rather than" injury.**
   - Female mice lose S3 identity genes on the same schedule (Dixon 2022).
   - Sham surgery causes the testosterone drop (Soljancic 2013).
   - AR is required for part of the injury response (Allison 2025).
   - The two processes co-occur and can be separated.

## Sex differences along the PT and androgen dependence

- **Ransick 2019** [full text], 2 mice per sex: "The highest disparity in gene activity between the sexes was in the
  S3 region". Cyp7b1 is male-restricted and Slc7a12 female-restricted in S3 (RNAscope). No hormone manipulation.
- **Xiong 2023, Dev Cell** [full text]: castration, testosterone, and Six2-, Sox2- and Alb-Ar knockouts, with
  snMultiome.
  - "After castration, 85.5% of the core M-biased genes (288 of 337) were down-regulated".
  - "On nephron-specific AR removal … 82.5% of core M-biased genes … were down-regulated".
  - "In male PT-S1, little effect was observed on either gene expression or chromatin status upon AR removal".
  - Slco1a1 is "entirely absent upon AR removal".
  - The scout counted AR-maintained genes in the public knockout tables (q < 0.05, |log2FC| ≥ 1):
    - S1 18, S2 200, S3 306;
    - AR-repressed: 25, 153 and 387.
  - Only about 30% of male-biased S3 genes are AR-dependent.
- **Chen L 2023, JASN** [full text], microdissection: sex-differential transcripts S1 79, S2 1,177, S3 766. No
  ARE/ERE over-representation in promoters. Ar mRNA is higher in females.
- **Chen S 2025, Nat Genet** [full text]: "Most sex-biased DEGs were found in S3, followed by S2"; the AR motif is
  enriched in males in all PT segments.
- **Harris 2021** [full text]: kidney-specific Ar knockout lowers PT cell height in PCT and cortical PST, not in the
  outer medulla.
- **KAP is atypical.** Its androgen-dependent compartment is cortical; S3 KAP is thyroid- and pituitary-driven
  (Meseguer & Catterall 1987 and 1990 [abstract only]; Teixidó 2006 [full text]). In the knockout tables Kap is S2
  −4.81, S1 −1.38 and S3 −1.15.

## Androgens around AKI

- **Soljancic 2013** [full text], rat:
  - "both sham rats and rats subjected to renal I/R had a 90% reduction in serum testosterone by 3 h after
    reperfusion";
  - recovery begins by 12 h in sham but not after I/R;
  - castration worsened IRI and testosterone at 3 h protected.
- **Patil 2016** [full text]: at 24 h I/R testosterone is 20 ng/dl against about 60% of baseline in sham. The dose
  effect is biphasic.
- **Robert 2011** [full text]: testosterone 2.4 → 0.7 nmol/L at day 1, with a reciprocal oestradiol surge.
- **Huang 2022** [full text]: renal AR mRNA falls in males at 24 h, and AR silencing worsens injury.
- **The direction is contested.**
  - Park 2004 (PMID 15358759) [abstract only]: castration protects, but "Antagonists of androgen or estrogen
    receptors do not affect the gender differences".
  - Kang 2014 [full text]: castration protects.
  - Müller 2002 [abstract only].

## Has anyone said this?

- **Kirita 2020** [full text]: male mice only. FR-PTC are defined by the loss of Slc5a12, Slc22a30 and Slc7a13 (an
  AR-dependent S3 gene) as differentiation markers. Its Dataset S3 SCENIC `Ar_extended` activity:
  - healthy S1 −0.46, S2 +0.80, S3 +1.97;
  - injured S1/2 −0.80, injured S3 −0.71, severe injured −0.89.
  - It ranks 7th of 134 regulons in healthy S3. Being AUCell on a 31-gene regulon, it is not independent evidence of
    lost signalling.
- **Gerhardt 2021** [full text]: "ectopic activation of Slc7a12, a female-specific S3-marker gene, in IRI replicates"
  (RNAscope). Hormones are not mentioned.
- **Gerhardt 2023** [full text]: Ar, Nr3c1 and Nr3c2 chromVAR activity rises in injured cells. These receptors share
  the response element, so a motif-based "AR activity" is confounded by glucocorticoids.
- **Allison 2025** [full text, bioRxiv]: "Castration or proximal tubule-specific deletion of Ar in males abolished
  Sox9 induction". PT Ar deletion reduced injury.
- **Ide 2022** [full text]: orchiectomy induces glutathione-metabolism genes in male kidneys.
- **Wu 2020** [full text], UUO: 64.8% of sex-differential PT genes change during UUO.
- **Dixon 2022** [full text], female mice: "many proximal tubule genes were rapidly downregulated after injury". Its
  S3-down set shares 10 genes with Gerhardt's male S3-down set, and neither is the AR program.
- **Viñas 2020** [abstract only]: its GEO DESeq2 contrasts (GSE144622) give IR against sham at 24 h in male and
  female sorted PT. Females respond much less (497 against 4,101 genes).

## Predictions carried into the protocol

1. **Feminisation.** AR-repressed genes rise while AR-maintained genes fall. If both fall, the change is not
   feminisation.
2. **AR dependence beyond the gradient and transport.** AR-maintained genes are 77% S1→S3-rising against 47% of
   others, so this is a real confound.
3. **Kirita timing.** An androgen withdrawal should open a gap that grows from 4 h to 12 h–2 d. The likely outcome is
   a maximal drop already at 4 h.
4. **A glucocorticoid-response control.**
