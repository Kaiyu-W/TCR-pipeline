# `TCR-pipeline`
This pipeline is designed to process `TCR-seq` data, with input data generating from `MiXCR` or `CellRanger` (the latter of which file format needs to transform by hand).

It needs `Java`, `Python`, `R`, each with extra packages `VDJTools`, `Numpy`/`Pandas`/`Matplotlib`/`Seaborn`, `DESeq2`/`ggplot2`. Please make sure Java/Python/Rscript works, and install the packages in advance (except VDJTools that has been embedding into this script), otherwise the pipeline will failed to work. Recommended to create `Conda` environment to manager and call the packages.

The way this script works is to unzip the temporary script into '/tmp/', and delete that at the end of the program.

More details please see the script. Please contact me if you need help. 
I'll complete the workflow in the future if I have more time as well as skills, since it lacks the code that helps to install the dependent packages.

Pipeline written by `Kaiyu Wang`, in SIBCB of CAS, Shanghai, in Sep 2021.
