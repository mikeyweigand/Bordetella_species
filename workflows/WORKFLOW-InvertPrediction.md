## Workflow for inversion prediction:  
For calculating a linear model of symmetric inversion breakpoints and then predicting all possible inversions between IS elements (or other repetitive sequences).  

1. Find all symmetric inversion pairs and calculate breakpoint distances.  
 + Follow __WORKFLOW-rearrangement.md__
1. Model and predict new boundaries
 + __inverts.predict.pl__  
1. For figure drawing see R code within:
 + __inverts.pred.R__
 + __circos-plots.R__     
