The HIPPIE network construction tool (HIPPIE_NC) takes as an input a flat file version of the HIPPIE database and and list of proteins or interactions around which a HIPPIE subnetwork is constructed. HIPPIE_NC is available from http://cbdm.mdc-berlin.de/tools/hippie/NC/HIPPIE_NC.jar .

Requirements and installation
=============================

It is written in java and only requires a recent version of the JDK installed. No further installation is necessary the tool can be run right away. The tool reads the interaction network HIPPIE from a flat file. The newest version can be found at: http://cbdm.mdc-berlin.de/tools/hippie/hippie_current.txt. As a default HIPPIE_NC expects the flat file hippie_current.txt to reside in the same directory as the script itself.


Quick Start
===========

The only obligatory parameter is the path to an input file listing either proteins or pairs of proteins. Assuming that the input file is named query.txt and is located in the same directory as the HIPPIE_NC script, the simplest command to run HIPPIE_NC is:
>java -jar HIPPIE_NC.jar -i=query.txt


Input data format
=================

Each line of the input file should contain either a single protein or a protein pair separted by either whitespaces or a tab. Unlike the web tool only uniprot ids and entrez gene ids are allowed. The type of identifier used throughout the input file has to be specified with the '-t' flag. As a default uniprot identifiers are assumed.

Options
=======
-h=HIPPIE (optional): path to HIPPIE (default assumption is that HIPPIE is in same directory and named hippie_current.txt)
-i=input: path to input file
-t=type (optional): either 'e' for entrez or 'u' for uniprot (default 'u').
-s=score: cutoff value for the interaction score between 0 and 1. All interactions with a score less then the indicated score are skipped (default 0.0).
-l=layer (optional): specifies if only interactions within the input set are considered (l = 0), between HIPPIE and the input set (l=1) or also neighbors of neighbors (l=2); maximum amount of layers is 3
-o=output (optional): if a file name is specified, the output is written to this file; otherwise to standard out (default standard out)


Example query
=============
>java -jar HIPPIE_NC.jar -i=query.txt -t=e -o=out.txt -s=0.5

For any questions, comments or suggestions, please send an email to martin.schaefer@crg.eu
