# m6Aboost
Author: You Zhou, Kathi Zarnack    

---

## Introduction
N6-methyladenosine (m6A) is the most abundant internal modification in mRNA. 
It impacts many different aspects of an mRNA's life, e.g. nuclear export, 
translation, stability, etc.   

m6A individual-nucleotide resolution UV crosslinking and immunoprecipitation 
(miCLIP) and the improved **miCLIP2** are m6A antibody-based methods that allow 
the transcriptome-wide mapping of m6A sites at a single-nucleotide resolution. 
In brief, UV crosslinking of the m6A antibody to 
the modified RNA leads to truncation of reverse transcription or C-to-T 
transitions in the case of readthrough. However, due to the limited specificity 
and high cross-reactivity of the m6A antibodies, the miCLIP data comprise a 
high background signal, which hampers the reliable identification of m6A sites 
from the data. 

For accurately detecting m6A sites, we implemented an AdaBoost-based machine 
learning model (**m6Aboost**) for classifying the miCLIP2 peaks into m6A sites 
and background signals. The model was trained on high-confidence 
m6A sites that were obtained by comparing wildtype and _Mettl3_ knockout mouse 
embryonic stem cells (mESC) lacking the major methyltransferase Mettl3. For 
classification, the m6Aboost model uses a series of features, including the 
experimental miCLIP2 signal (truncation events and C-to-T transitions) as well 
as the transcript region (5'UTR, CDS, 3'UTR) and the nucleotide sequence in a 
21-nt window around the miCLIP2 peak.

The package [m6Aboost](http://bioconductor.org/packages/m6Aboost) includes the 
trained model and the functionalities to prepare the data, extract the 
required features and predict the m6A sites.

---

## Citing m6Aboost

> Nadine Körtel, Cornelia Rücklé, You Zhou, Anke Busch, Peter Hoch-Kraft, F X
> Reymond Sutandy, Jacob Haase, Mihika Pradhan, Michael Musheev, Dirk Ostareck,
> Antje Ostareck-Lederer, Christoph Dieterich, Stefan Hüttelmaier, Christof
> Niehrs, Oliver Rausch, Dan Dominissini, Julian König, Kathi Zarnack, Deep
> and accurate detection of m6A RNA modifications using miCLIP2 and m6Aboost
> machine learning, Nucleic Acids Research, Volume 49, Issue 16, 20 September
> 2021, Page e92, https://doi.org/10.1093/nar/gkab485

---

## How to use it
Documentation (vignette and user manual) is available at the **m6Aboost's** 
Bioconductor landing page at http://bioconductor.org/packages/m6Aboost.

---


