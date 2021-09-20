# TCR-pipeline
This pipeline is designed to process `TCR-seq` data, with input data generating from `MiXCR` or `CellRanger` (in which file format needs to transform by hand).

It needs `Java`, `Python`, `R`, each with extra packages `VDJTools`, `Numpy`/`Pandas`/`Matplotlib`/`Seaborn`, `DESeq2`/`ggplot2`. Please make sure Java/Python/Rscript works, and install the packages in advance (except VDJTools that has been embedding into this script). Recommended to create `Conda` environment to manager and call the packages.

The way this script works is to unzip the temporary script into '/tmp/', and delete that at the end of the program.

More details please see the script. Please contact me if you need help.

Pipeline written by `Kaiyu Wang`, in SIBCB of CAS, Shanghai, in Sep 2021.
