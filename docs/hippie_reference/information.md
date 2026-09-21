## HIPPIE Howto

- [Generation of reliable and meaningful PPI networks with HIPPIE.](#welcome)

- [Confidence scoring of experimentally measured interactions.](#scoring)

- [Query interactions.](#query)

- [Network construction.](#network)
   

      
  - [Input.](#input)
      
  - [Output type.](#output)
      
  - [Interaction layers.](#int_layers)
      
  - [Min. number of PPIs to query set.](#common_int)
      
  - [Score filter.](#score)
      
  - [Interaction type filter.](#type)
      
  - [Tissue filter.](#tissue)
      
  - [Functional filter.](#function)
      
  - [Directed edges.](#direction)
      
  - [Inhibitory or activating effect of an interaction.](#effect)
   

- [Exporting query results.](#export)

- [Disease and Gene Ontology enrichment of query results.](#enrichment)

- [Browsing the HIPPIE proteome.](#browsing)

- [Cytoscape visualization of HIPPIE output.](#cytoscape)

- [HIPPIE network construction script.](#script)

- [Annotation of a list of newly measured interactions.](#annotation)

- [Download the repository.](#download)

- [Example of a HIPPIE workflow.](#example)

- [Sources.](#sources)

- [HIPPIE API.](#api)

- [General technical remarks.](#remarks)

- [Related manuscripts (how to cite).](#cite)

- [Contact.](#contact)

#### Generation of reliable and meaningful PPI networks with HIPPIE
Welcome to HIPPIE (the Human Integrated Protein Protein Interaction Reference)! We provide confidence scored and functionally annotated human protein-protein interactions (PPIs). This page contains details on the information stored in HIPPIE and the multiple ways to access it. First time users might want to start with the [walk-through example](#example) to get an idea what types of analyses can be performed with HIPPIE.

**[Back to top](#top)

#### Confidence scoring of experimentally measured interactions
A core component of HIPPIE is the confidence scoring of interactions based on the amount and reliability of evidence supporting them.
This score is calculated as a weighted sum of the number of studies in which an interaction was detected, the number and quality of experimental techniques
used to measure an interaction and the number of non-human organisms in which an interaction was reproduced.

The parameters of this scoring scheme were jointly optimized by a group of experts and a computer algorithm: we first assigned quality scores to each experimental technique that measures protein-protein interactions.
This experimental quality score is supposed to reflect the reliability and the error rate of the techniques. The list of experimental quality scores
can be found in the [download](download.php) section. We are aware that, to a high degree, the estimation of experimental quality is led by the individual perception and experience of our expert team and can
only reflect our subjective beliefs. If you feel that a technique was fundamentally misjudged, please give us your feedback.

Also, there is the possibility for users of HIPPIE to modify the scores and to apply our rescoring tool to get an individually scored version of HIPPIE. Please, find a stand-alone version of the tool and instructions on
how to run it in the [download](download.php) section.

**[Back to top](#top)

#### Query interactions
The [query mode](index.php) allows for querying HIPPIE with a [UniProt](http://www.uniprot.org/) identifier (id or accession), [gene symbol](http://www.ncbi.nlm.nih.gov/gene) or [Entrez](http://www.ncbi.nlm.nih.gov/gene) gene id.
All literature interactions in HIPPIE that match the query protein are displayed.

On the result page, UniProt id, Entrez gene id and gene symbol of the interactors of the query protein are indicated.
The confidence score of the interaction is given in the last column. Clicking on it opens an evidence page where the source databases of the interaction are listed, together with
the publications which mention the interaction, the experimental systems used to detect it and the species in which it is conserved.

**[Back to top](#top)

#### Network construction
In the [network construction mode](network.php), several proteins (or interactions) can be queried at the same time ("batch mode") and a network is constructed from the input. The resulting network
can be either visualized graphically, listed in a tabular format in the browser or downloaded as a file in several formats.

##### Input
In the network construction mode HIPPIE can be either queried with a list of proteins or a list of interactions (a maximum of 100 lines is allowed, that can be either pasted into the query field or uploaded as a file).
Each line must contain a single protein or a protein pair (separated by a tab or whitespaces). Allowed identifier types are: [UniProt](http://www.uniprot.org/) identifiers (id or accession),
[gene symbols](http://www.ncbi.nlm.nih.gov/gene) or [Entrez](http://www.ncbi.nlm.nih.gov/gene) gene ids.

The main goal of this mode is to allow users to construct networks around proteins of interest or around experimentally measured protein interactions.
For that purpose also interacting protein pairs can be used as seeds for the network construction and both the input interactions and interactions from HIPPIE
will be displayed in the resulting network. Proteins and interactions submitted by the user are highlighted in blue in our graphical output and, in case they are not in HIPPIE, they are highlighted in red.
Several filter criteria allow to modify the resulting network (see below). Note, however, that these filters are only applied to interactions added from the HIPPIE repository and do not affect
the input proteins or interactions, which will be displayed regardless of whether they pass the filter criteria or not.

**[Back to top](#top)

##### Output type
The user can choose between different output types.

    
- show in browser - text: output the table of interactions in the browser
    
- show in browser - visualization: use Cytoscape.js to visualize the constructed subnetwork
    
- HIPPIE tab file format : the interactions are written to a simple text file
    
- PSI-MI TAB 2.5 output format : interactions are written to the standardized MITAB format

**[Back to top](#top)

##### Interaction layers
The size of the network can be modulated by choosing the interaction layer size:

    
- 0 - only interactions within the input set are considered
    
- 1 - all interactions in which members of the input set participate are listed

Please note that we had to limit the amount of layers that can be queried and the size of the input sets due to the computing time that rapidly increases for large input sets or layer-2 queries.
However, we offer a script for [download](download.php) that implements the algorithm behind the HIPPIE network construction tool and can be run locally.

**[Back to top](#top)

##### Min. number of PPIs to query set
This output parameter is useful when the user wants to know whether the query set has one or more common interactors. The default value of 1 simply queries for proteins interacting with the query set.
However, a value of 2, for example, focuses on proteins that interact with at least 2 members of the query set.

Note that only integer values greater or equal to 1 are valid. Values below 0 are reset to 1 and fractions are rounded.

**[Back to top](#top)

##### Score filter
A threshold on the HIPPIE confidence score can be chosen. The user can either specify a custom value between 0 and 1 or choose a predefined confidence level:
medium confidence (0.63 - second quartile of the HIPPIE score distribution) or high confidence (0.73 - third quartile). If both custom and predefined thresholds are selected, the higher one is applied.

Note, that the confidence filter is not applied to interactions uploaded by the user.

**[Back to top](#top)

##### Interaction type filter
HIPPIE offers an interaction type filter to distinguish between binary and complex interactions. For the interaction type filter we rely on the annotation of interactions from MINT, BioGRID and IntAct.
We offer to select interactions annotated with the PSI-MI categories association (MI:0914), physical association (MI:0915), direct interaction (MI:0407) and colocalization (MI:0403).

**[Back to top](#top)

##### Tissue filter
Interactions can be filtered for tissue specificity. Several tissue filters can be selected at the same time by ticking more than one tissue from the menu.
For an interaction to pass this filter, both interactors should be expressed in at least one of the selected tissues. 

It is also possible to upload user defined filter sets (e.g. results from microarray experiments), which are used instead (if no tissues are selected) or in addition (if there are selected tissues) to
the predefined tissue filter sets.

The tissue expression RNA-Seq data was taken from [GTEx](http://gtexportal.org/). A gene was considered to be expressed in a given tissue if it showed an RPKM >= 1.

##### Functional filter
Interactions can be filtered for annotation of [GO](http://www.geneontology.org/) and [MeSH](http://www.nlm.nih.gov/mesh/) terms. Interactions are associated
with the lowest common ancestor within the hierarchy of terms associated with the interacting proteins. Selection of a more generic term causes interactions annotated with a more specific term to pass the filter. 

**[Back to top](#top)

##### Directed edges
HIPPIE infers edge directionality by calculating the shortest paths that connect sources with sinks to infer information flow in signaling networks. The source and sink sets can be specified by pasting gene lists into the respective input fields. Receptors and transcription factors (according to GO annotations) can be selected as default source and sink sets. The shortest path computation is done using the python package NetworkX (https://networkx.github.io/). The user can choose between two ways of shortest path computation: 1. All shortest paths consisting of the same number of edges (not considering edge weights) are highlighted ('show unweighted shortest paths'). 2. If there are several shortest paths consisting of the same number of edges, only the one with the highest cumulative confidence score is shown ('show confidence-weighted shortest paths').

Alternatively, for known pathways, KEGG's directionality assignment can be displayed. Edge directionality is indicated in the visual output mode with arrows. Sources and sinks are color-coded. In the tabular output mode (browser or download) an additional column indicates the direction. 

**[Back to top](#top)

##### Inhibitory or activating effect of an interaction
We predicted the inhibitory and activating effect of interactions based on phenotypic similarity upon gene knockdown ([Suratanee et al.](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003814)). We could predict an effect for 8316 interactions in HIPPIE. This includes several newly included interactions from databases that are not limited to experimentally verified interactions or do not provide experimental evidence. These interactions have been assigned a confidence score of 0. To include these interactions, a custom score threshold of 0 must be set (see [Score filter](information.php#score_filter)). Additionally, we integrated effect assignments from KEGG. The effect information is visually encoded in the graph visualization output mode by edges terminating in either arrows (activation) or bars (inhibition). If the effect visualization is chosen in combination with edge directionality, a diamond-shaped arrow indicates that an edge has a direction but not an effect. Edges with two arrows have an effect but not a direction. In addition to the visual encoding, the interaction effects can be shown in the tabular output. 

**[Back to top](#top)

#### Exporting query results

The browser-based results of the above-mentioned [protein](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/index.php) and [network](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/network.php)
queries can be copied or downloaded as tab-separated text files (TSV format).
When the network visualization output is chosen in the network query, the resulting Cytoscape.js representation of the subnetwork can be saved as a PNG or JPEG graphic, or as a JSON object.

**[Back to top](#top)

#### Disease and Gene Ontology enrichment of query results

It is possible to perform a Disease Enrichment Analysis with the proteins resulting from a [protein](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/index.php) or
[network](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/network.php) query. The set of diseases that are significantly associated to the partners of a protein or the members of a subnetwork is shown,
respectively. This analysis is performed with the help of the text mining tool [Gene set to diseases](http://cbdm-01.zdv.uni-mainz.de/~jfontain/cms/?page_id=592), which computes a disease enrichment
analysis on gene sets, using biomedical literature data ([click here](http://cbdm-01.zdv.uni-mainz.de/~jfontain/cms/?page_id=592) for more details).

It is also possible to perform a Gene Ontology Enrichment Analysis of the proteins resulting from a [protein](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/index.php) or
 [network](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/network.php) query. The analysis is performed via [PANTHER](http://pantherdb.org/), the preferred tool of the
 [Gene Ontology Consortium](http://geneontology.org/) ([click here](http://nar.oxfordjournals.org/content/33/suppl_1/D284.full) for more details).

These analyses are restricted to 1000 proteins.

**[Back to top](#top)

#### Browsing the HIPPIE proteome.

The [browsing mode](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/browse.php) in HIPPIE allows for the exploration of the entire set of proteins in the database and presents summary
statistics about the network members. UniProt, Entrez and Symbol IDs, node degree (i.e. the number of interaction partners), average score of adjacent edges and a link to query interaction partners are shown
for each protein in a big table. Note that loading the entire HIPPIE proteome may take a few seconds, depending on the type of Internet connection being used.

**[Back to top](#top)

#### Cytoscape visualization of HIPPIE output
HIPPIE offers to display any constructed subnetwork with the light-weight network visualization tool [Cytoscape.js](http://js.cytoscape.org/).
However, for more detailed analyses, we recommend using HIPPIE's export facilities and using the desktop version
of [Cytoscape](http://www.cytoscape.org/). The HIPPIE tab format is specifically designed in a way that
important features of interactions are listed in separated columns and can be easily imported into Cytoscape as edge features.
For that purpose choose 'File' -> 'Import' -> 'Network from table' in Cytoscape. Select a file exported from HIPPIE and activate the columns which should be imported.

In our example the 'score' and the 'uploaded interaction' columns have been selected and are imported as edge features. In the resulting Cytoscape network (see below) the score attribute is size encoded:
high confidence interactions are drawn with a thicker line as compared to those with fewer experimental evidence. The information whether an interaction is from HIPPIE or uploaded by the user is color encoded:
HIPPIE interactions are painted in blue while user specified interactions are rendered in red.

Alternatively, a visual representation of the network can be generated using HIPPIE and Cytoscape.js and then exported for further analyses. The available JSON format can be imported into Cytoscape
(requires Cytoscape version 3.1.0 or newer).

**[Back to top](#top)

#### HIPPIE network construction script
For large input sets or queries involving indirect neighbors (layer-2 or -3 interactions), we developed a script that can be downloaded and run locally.
It is implemented in Java, making it platform-independent. Instructions on how to use it can be downloaded together with the tool from the [download](download.php) page.

**[Back to top](#top)

#### Annotation of a list of newly measured interactions
HIPPIE aims at serving experimentalists conducting large scale screens for protein-protein interactions.
Hence, it offers a service in which a result table from such a screen is [annotated](annotation.php) with information
on which interactions have been observed before and how reliable they are under HIPPIE's scoring scheme.
The input format is very flexible, allowing the user to select from different ID types (Uniprot and Entrez),
indicate the columns where the interactors are listed and specify parameters for the input format
(which characters separate columns and whether there is a preceding header line).

**[Back to top](#top)

#### Download the repository
HIPPIE is provided in two different download formats: PSI-MITAB 2.5 and our own tab-separated flat file format.
Besides offering the data in compliance with the specification of the [Proteomics Standards Initiative (PSI)](http://www.psidev.info/),
we decided to additionally offer the data in a format described below which makes it slightly easier to either
compare or extend experimental datasets with HIPPIE, by resolving ambiguities and synonymous ids into several lines: 

The columns indicate (1) UniProt identifier and (2) Entrez Gene identifier of the first protein partner,
(3) UniProt identifier and (4) Entrez Gene identifier of the second protein partner, (5) score and (6) a description summarizing the origin of the evidence for the interaction.
If one gene maps to several proteins each combination of proteins is listed in a separate line.

**[Back to top](#top)

#### Example of a HIPPIE workflow
To demonstrate HIPPIE's functionality, we give a simple example of a typical query that makes use of expression data and shortest paths computations.

We start from the NETWORK QUERY tab. The user can enter a list of proteins or protein pairs (separated by a space or tab, which will then be interpreted as interacting proteins) into the up-most text field
(encircled below). In this example, we want to explore HIPPIE's functionality to construct context-specific signaling networks. We therefore query HIPPIE with the kinases BRAF, MEK1 (MAP2K1) and ERK1 (MAPK3),
which are members of the Mitogen-activated protein kinase (MAPK) signaling cascade and activate each other in the stated order. 

Additionally, we select two filters: (A) high confidence interactions to make sure that only experimentally reliable interactions are shown. (B) Colon tissue to exclude all PPIs formed by proteins not expressed
in colon. The activity of the MAPK pathways is altered in many colon cancers. 9% of all colon cancer patients carry BRAF mutations (Lawrence et al. Nature. 2014).
(C) We enable shortest path computation from BRAF to transcription factors

Displaying shortest paths between BRAF ("source", in blue) and transcription factors ("sinks") correctly reproduces the chain of signaling events (BRAF activates MEK1, which activates ERK1 in turn).
All terminal nodes in pink (ELK1, MYC, JUN, TP53, SREBF1/2) are known substrates of ERK1 (Yoon & Seger. Growth Factors. 2006)Γüá.

The resulting network can be further analysed. One option is to identify diseases overrepresented among the proteins in the network. To this end, the proteins in the network will be automatically transferred to
the tool "Gene set to Diseases" (Fontaine & Andrade-Navarro. Under review) when the link in the circle is clicked.

In the case of the MAPK example, the disease enrichment analysis reveals that the subnetwork is strongly enriched in ΓÇ£neoplastic cell transformationΓÇ¥ (q < 1e-9). In agreement with filtering for colon expression,
the most strongly represented single tumor type is 'colonic neoplasms' (q < 1e-5).

**[Back to top](#top)

#### Sources
HIPPIE integrates interaction data from 10 source databases and 11 studies (that have not been fully covered by the other databases yet). The chart below visualizes the contributions
of the main sources. Only studies and databases contributing with at least 200 unique (not found in any other dataset) interactions are shown.
 
 [Bell09](http://www.ncbi.nlm.nih.gov/pubmed/19293945) is a high-throughput study.
 The bar plot was created using [Google chart tools](http://code.google.com/apis/chart/).

**[Back to top](#top)

#### HIPPIE API
HIPPIE users can query the resource via our REST web service. This means that HIPPIE can be easily integrated in any Bioinformatics pipeline using the following template:

http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/queryHIPPIE.php?**proteins=**xxx,xxx;xxx|xxx&**layers=**xxx&**conf_thres=**xxx&**out_type=**xxx

   
- **proteins** = One or more proteins of interest separated by ",", ";" or "|" (*mandatory*).
   
- **layers** = 0 to query interactions within the input set or 1 to query interactions between the input set and HIPPIE (*optional, default = 1*).
   
- **conf_thres** = Only protein interactions with confidence scores above this threshold, which can range between 0 an 1, are considered (*optional, default = 0*).
   
- **out_type** = The query output format. *browser* shows the list of interactions in a table in HIPPIE, *viz* shows a network visualization, *mitab* generates a MITAB file and *conc_file* generates a simple tab-separated text file (*optional, default = conc_file*).

**[Back to top](#top)

#### General technical remarks
HIPPIE works best with a recent version of Mozilla Firefox (≥3.6.). We observed problems with Microsoft Internet Explorer. While we are fixing these problems, we recommend you to visit HIPPIE with a different web browser.

**[Back to top](#top)

#### Related manuscripts

- Schaefer MH, Fontaine J-F, Vinayagam A, Porras P, Wanker EE, et al. (2012) *HIPPIE: Integrating Protein Interaction Networks with Experiment Based Quality Scores*. PLoS ONE, 7(2): e31826 [[link](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0031826)]

- Schaefer MH, Lopes TJS, Mah N, Shoemaker JE, Matsuoka Y, et al. (2013) *Adding Protein Context to the Human Protein-Protein Interaction Network to Reveal Meaningful Interactions*. PLoS Computational Biology, 9(1): e1002860 [[link](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002860)]

- Suratanee A, Schaefer MH, Betts M, Soons Z, Mannsperger H, et al. (2014) *Characterizing Protein Interactions Employing a Genome-Wide siRNA Cellular Phenotyping Screen*. PLoS Computational Biology, 10.9: e1003814 [[link](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003814)]

- Alanis-Lobato G, Andrade-Navarro MA, & Schaefer MH. (2016). *HIPPIE v2.0: enhancing meaningfulness and reliability of protein-protein interaction networks*. Nucleic Acids Research, gkw985
[[link](http://nar.oxfordjournals.org/content/early/2016/10/28/nar.gkw985)]

**[Back to top](#top)

#### Contact

Please send questions or comments on HIPPIE to: [martin.schaefer@ieo.it](mailto:martin.schaefer@crg.eu)
HIPPIE is hosted by the [CBDM group](http://cbdm.uni-mainz.de/) at the JGU in Mainz.

**[Back to top](#top)