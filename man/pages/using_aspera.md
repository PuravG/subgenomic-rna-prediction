# Using Aspera

The IBM Aspera connect tool can increase download speeds for large genomic datasets from SRA. This is free but proprietary software from IBM, and is not available via conda or pip.


To use the tool, download the appropriate installer for your machine [directly from IBM](https://www.ibm.com/support/fixcentral/swg/selectFixes?parent=ibm~Other%20software&product=ibm/Other+software/IBM+Aspera+Connect&release=All&platform=All&function=all). **Please download a release no later than 4.1.3.93** All releases from 4.2.0 onwards do not provide an openssh key, which is a requirement for the enaDataGet tool we use to get data from RNA.


You will then need to provide two file paths in your config.ini file.

`aspera_bin` is the path to the folder containing the executable `asperaconnect`, in Linux by default this is something like

`/home/yourname/connect/bin`

`aspera_openssh` is the path to the file `asperaweb_id_dsa.openssh`, something like

`/home/yourname/connect/etc/asperaweb_id_ds.openssh`