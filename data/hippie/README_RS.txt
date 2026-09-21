The HIPPIE rescoring tool (HIPPIE_RS) takes as an input a flat file version of the HIPPIE database and a list of scores assigned to experimental techniques. HIPPIE_RS is available from http://cbdm.mdc-berlin.de/tools/hippie/RS/HIPPIE_RS.jar .

Requirements and installation
=============================

The script is written in java and only requires a recent version of the JDK installed. No further installation is necessary, the tool can be run right away. It reads the interaction network HIPPIE from a flat file. The newest version of HIPPIE can be found at: http://cbdm.mdc-berlin.de/tools/hippie/hippie_current.txt. As a default HIPPIE_RS expects the flat file hippie_current.txt to reside in the same directory as the script itself.


Quick Start
===========

When HIPPIE and the list of experimental scores are in the working directory, the tool can be run without passing any further parameters. 
>java -jar HIPPIE_RS.jar
A rescoring run putting a higher emphasis on the orthology subscore but leaving the saturation parameters as they are (for details on the scoring formula please check out the homepage or the HIPPIE paper) could look like: 
>java -jar HIPPIE_RS.jar -ws=0.6 -wo=0.1 -wt=0.3

Options
=======
-h=HIPPIE: path to HIPPIE (default assumption is that HIPPIE is in same directory and named hippie_current.txt)
-e=input: experimental score file (default assumption is that this is in same directory and named experimental_scores.tsv)
-o=output: if a file name is specified, the output is written to this file; otherwise to standard out (default standard out)
-as: parameter controling how fast the study subscore is saturating (between 0 and infinity), default value (as used in the paper) is 2.3
-ao: parameter controling how fast the orthology subscore is saturating (between 0 and infinity), default value (as used in the paper) is 1.6
-at: parameter controling how fast the technique subscore is saturating (between 0 and infinity), default value (as used in the paper) is 0.2
-ws: parameter controling the contribution of the study subscore to the overall score (between 0 and 1), default value (as used in the paper) is 0.6
-wo: parameter controling the contribution of the orthology subscore to the overall score (between 0 and 1), default value (as used in the paper) is 0.1
-wt: parameter controling the contribution of the technique subscore to the overall score (between 0 and 1), default value (as used in the paper) is 0.3 (note that the three weights ws, wo and wt have to sum up to 1).


For any questions, comments or suggestions, please send an email to martin.schaefer@crg.eu
