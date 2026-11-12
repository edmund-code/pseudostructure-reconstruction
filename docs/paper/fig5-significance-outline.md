# Fig. 5 Results: outline with significance per result

Working outline (2026-11-11). Each result: what we measured, why it matters clinically, what follows.
Sources were read in full text; quotes are in `citation_ledger.md`. "Speculate" marks anything beyond our data.

---

## R1. SGLT2 occupies a longer stretch of the human PT

- **Our data.** *Slc5a2* above half-maximum through 0.38–0.39 of the mouse PT, 0.66–0.71 of the human PT.
  *Slc2a2* (GLUT2) ends at the same place in both (0.55–0.58). Both swaps < 1/4 of the contrast.
- **Significance.** SGLT2 is the gliflozin target; gliflozins are approved for type 2 diabetes, heart failure and
  chronic kidney disease [PMC11462946]. SGLT2 carries ~97% of filtered glucose; under inhibition, SGLT1 downstream
  recovers part of it, leaving 50–60% excreted in mice [Rieg 2014]. Rat modelling: inhibition moves Na+ reabsorption
  downstream, S3 oxygen consumption +102% [Layton 2015].
- **Implication.** Both depend on how much PT lies past SGLT2. Longer human SGLT2 domain -> shorter downstream
  stretch -> mouse estimates of SGLT1 compensation and late-PT workload may not carry over; also relevant to dual
  SGLT1/2 inhibition (sotagliflozin).
- **Status.** Written and committed (e7868f6).

---

## R2. Secretory transporters sit at different positions

- **Our data.** OAT1 (*Slc22a6*) and OAT3 (*Slc22a8*) peak at the human PT start (0.06–0.16), early-to-mid in mouse
  (0.24–0.47). MRP2 (*Abcc2*) starts at the human PT start, past the mouse midpoint. OCT2 (*Slc22a2*) starts later in
  human (0.51–0.57 vs 0.32–0.36). URAT1 (*Slc22a12*) peaks at 0.84 in human, 0.54 in mouse. All swaps < 1/4.
- **Significance (sources checked, text still to write).**
  - OAT1/OAT3: uptake of anionic drug classes; regulatory drug-drug interaction transporters; coupled to apical
    MRP2/MRP4 efflux [Nigam 2015, PMC4281586]. Replaces the placeholder PMC3675737.
  - OAT1 and the antiviral tenofovir [Kohler 2011, PMC3103636].
  - OCT2 and the chemotherapy agent cisplatin; human *SLC22A2* variant [Filipski 2009, PMC2746866].
  - URAT1: urate reabsorption, target of uricosurics (benzbromarone, lesinurad) [PMC11811179].
- **Implication (to write).** Uptake and efflux co-localize at the human PT start but are offset in mouse; the human
  OCT2 window is shorter and later; the uricosuric site of action sits later in human. Mouse-based maps of where
  these drugs act should be checked by position before transfer.
- **Status.** Significance paragraph NOT written (three classifier stops). Sources verified.

---

## R3. Switch points move later along the human PT

- **Our data.** 53 pairs with a switch in all four sections; 13 (25%) shifted, 0 under swaps; 8 later in human.
  SGLT2 -> SGLT5 +0.21; SLC22A23 -> URAT1 +0.25; SLC22A23 -> OCT2 +0.24; y+LAT1 -> AGT1 +0.24; B0AT1 -> B0AT3 +0.18.
  7/7 checkable pairs reproduce in Lake 2023, including SGLT2 -> SGLT1.
- **Significance.** A later switch = the early partner carries more of the length in human. Matters when the early
  partner is a drug target. B0AT1 (*Slc6a19*) is a phenylketonuria target: *Slc6a19*-null mice excrete >100-fold more
  phenylalanine and the mouse model improves [Belanger 2018, PMC6124451].
- **Implication (to write).** An inhibitor of an early partner may cover a larger share of reabsorption in human than
  the mouse predicts (B0AT1); the later SGLT2 -> SGLT1 switch agrees with R1.
- **Status.** To write; source verified.

---

## R4. Segment-level clustering misses most of the shift

- **Our data.** Straight lines through S1/S2/S3 cluster means recover +0.07 of +0.25 (URAT1), +0.06 of +0.24 (OCT2),
  +0.02 of +0.21 (SGLT5); no switch for 3 pairs; one monocarboxylate pair moves the wrong way. Post hoc.
- **Significance.** Several switches move inside S2, which three means cannot resolve.
- **Implication (to write).** Drug-target maps built from segment-resolved atlases would place these handovers at
  similar positions in both species and miss the differences in R1-R3. One caveat already in the text: three cluster
  means are a simple baseline; finer clustering may recover more.
- **Status.** Result and caveat written; the "what this costs a reader" sentence still to add.

---

## R5. Mouse-expanded gene families carry extra zonation

- **Our data.** Mouse: 9.9% (1:many) and 14.2% (many:many) strongly zonated vs 3.8% (1:1); MH OR 2.24, P = 8e-7.
  No excess where the human counterpart expanded (1.7%). Human: same ordering, OR 1.57, P = 0.17.
- **Significance.** Families named in the text (*Cyp4a*, *Cyp2d*, *Ugt2b*, *Slco1a*) are drug-metabolizing or
  drug-transporting families [TO VERIFY - no source read yet].
- **Implication (to write).** Positional drug handling measured in the mouse PT may partly sit in gene copies with no
  single human counterpart; a one-to-one comparison drops them silently.
- **Status.** To write; needs one source for the "these families handle drugs" claim.

---

## Close

- COAST compares positional expression across species on one continuous coordinate.
- In the PT: glucose, secretory and urate transporters occupy different positions; switch points shift; segment-level
  clustering misses most of it.
- Limits already in the text: one donor, two cortical sections; transcript position != protein position; staining in
  more donors would confirm.

---

## Open

- R2 significance paragraph blocked three times by the response classifier; user to draft, Claude to check against
  sources and place.
- R5 needs a source for the drug-family claim.
- Replace placeholder PMC3675737 with Nigam 2015 (PMC4281586).
- Two different "Extended Data Fig. X" placeholders (steepness panel; *Abcc2* and coverage).
