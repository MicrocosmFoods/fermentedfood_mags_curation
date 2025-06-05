# Uploading Genomes and Sample Set Information to KBase

This document describes setting up a KBase narrative with a set of curated genomes and metadata for downstream users to explore and analyze. You can sign-up for a KBase account [here](https://narrative.kbase.us/legacy/signup) with any email account. You are not required to be affiliated with an academic institution or government organization to use KBase, and you can reside outside of the United States. You can find KBase documentation [here](https://docs.kbase.us/).

## About KBase
[KBase](https://www.kbase.us/) is a community-driven platform for facilitating open science research in systems biology. KBase allows you to run bioinformatics tools on large datasets using freely available Department of Energy high-perofrmance computing resources, allowing for open-sharing of research outputs and collaborative work. 

We have deposited our set of ~4000 "strain-representative" (clustered at 99% average nucleotide identity) genomes into KBase. This document describes how that process was done, and how users can also use these steps to create their own custom narratives containing either their own genomes, subsets of the strain-representative collection, or additional genomes from our database that are available on [Zenodo](TODO).

This documentation follows the bulk import KBase documentation [here](https://docs.kbase.us/data/upload-download-guide/uploads) with added steps for implementing a python script to help facilitate uploading genomes in batches.

## Create a Narrative
The main working unit on KBase is the narrative. A narrative is where you can create shareable, reproducible workflows containing certaint tools chained together, analyses and results, and/or datasets. 

To create a narrative, click on the `+ New Narrative` button in the top-right corner of the KBase interface. 
![Create Narrative](../img/create-narrative.png)

## Create the Bulk Import Template for your Batches
If you have a large set of genomes to upload (say about more than 1000), it is recommended to prepare them in batches. Since we had ~4000 genomes to upload, we split them into batches of around 500 genomes each. 

If you don't need to create the bulk import template (as in you can just use the python script out of the box) skip to the next step.

First to show how the bulk import works, you can upload a single genome to KBase. On the left-hand top pane of KBase, click `Add Data`. Then on the pop-out click `Import`. 

<div style="display: flex; gap: 10px;">
  <img src="../img/add-data.png" alt="Add Data" width="300" height="200">
  <img src="../img/import.png" alt="Import" width="300" height="200">
</div>


You can either upload data from a URL, Globus, or from your local computer. In our case, we uploaded FASTA files for genome assemblies from a local computer. Then select the import type as "FASTA Assembly". Then make sure to select that file and click "Import". 
![Import Type](../img/import-type.png)

Once you click import, a popout will appear in your narrative called `Import from Staging Area`. This will show the file(s) you selected from the staging area. If you choose to prep your batches of genomes using the python script below, you should only have to do this next step once. 

Fill out the two fields selecting the options that match your genome assembly. For our purposes, we are either uploading metagenome-assembled genomes (MAGs) or isolate genomes. For this example we selected metagenome-assembled genomes, and a minimum contig length of 500 bp. For our purposes we just selected a minimum contig length of 500 bp for every input assembly, but if you want to programatically enter this information that's entirely possible to implement with our python script documented below. 
![Create Template](../img/create-template.png)

Then click `CREATE IMPORT TEMPLATE` at the top right of the module. Then go back to "My Data" and your staging area. This template file will be in a folder called `bulk_import_templates`. Open the folder, and download the file. If the file doesn't automatically appear click "Refresh" a couple of times. 

![Bulk Template](../img/template-file.png)

Your bulk template CSV is going to look like this: 
```
Data type: assembly; Columns: 4; Version: 1			
staging_file_subdir_path	type	min_contig_length	assembly_name
FASTA file path	Assembly type	Minimum contig length	Assembly object name
YouL_2022__SRR13308542__bin.20.fa	Metagenome	500	YouL_2022__SRR13308542__bin.20
```

The python script below handles both setting up batches of genome files for upload and creating template specification import files that adhere to this structure.


## Prepare Batches of Genomes
The python script in `../scripts/batch_fasta_files` handles both programatically creating batches of genomes into subdirectories from a single input directory of genomes and creating input specification template files that correspond to those batches. 
```
usage: batch_fasta_files.py [-h] [--batch-size BATCH_SIZE]
                            [--output-dir OUTPUT_DIR]
                            [--samplesheet-dir SAMPLESHEET_DIR]
                            input_dir
batch_fasta_files.py: error: the following arguments are required: input_dir
```

For example, we used this script to split our input directory of ~4000 "strain-resolved" genomes into batches of 500 using this script. **Note that we configured this script so that assembly files starting with GCA/GCF are designated as "Draft isolates" and all other assembly files are designated as Metagenome-assembled genomes".** This might not be the case for your input files, so adjust the script as needed. 

We created the batches launching the script as follows: 
```
python3 scripts/batch_fasta_files.py 99id_genomes --batch-size 500 --output-dir 99id_batches --samplesheet-dir 99id_batch_samplesheets
```

In general it's probably not recommended to upload more than 500 at at time, but I'm not sure if this is a recommendation or hard and fast requirement. 

## Upload Genomes to KBase 
Now you are ready to upload each batch of genomes to KBase. When you navigate back to the "Data" pane and upload from your computer, select the first folder containing your batch of genomes. To select all files on a Mac press Cmnd + A. Then click Upload. **Do not refresh the page, the upload will be interrupted.**

Once all the genomes are uploaded, then upload the corresponding bulk template CSV. The script names the folders and the CSV files the same, such as `batch_1/` folder of genomes and `batch_1.csv` for the corresponding template CSV. Upload the CSV. 

If the batch CSV doesn't appear automatically, refresh the data in the staging area or refresh the Narrative (only if all your files have finished uploading to the staging area!!). 

Select the import CSV and select the file type as "Import Specification" and click "Import Selected". **Only import the CSV specification file, not any of the genomes!!**

![Bulk Import Specification](../img/bulk-import-specification.png)

A popout module "Import from Staging Area" will appear. Now all of your FASTA files in that batch appear without you having to manually click them all and import! (Note that the below shows files that have already been uploaded and shows a warning for this). If everything looks good, click "Run"
![Import Files](../img/import-files.png)

Now you will have to do this process for however many batches you have (in our case around 10). As your data uploads, it will start to appear in the left-hand side "Data" panel. 

## Upload Sample Set Information to KBase 

The associated metadata will have to be modified slightly to upload to KBase as a Sample Set. The script `sctipts/modify-metadata.py` handles making these small modifications, namely changing the first column to be "name" and changing long lists of sample or run accessions to lists of ranges if the samples/runs are lists that are in sequential order. Upload this to KBase in the same fashion as the genomes in the staging area, and set the file type to "Sample Set". This will need to be done so the metadata is only for the filtered set of genomes because the upload will crash with ~13,000 rows. Supposedly it will crash when you are over 10,000 rows, so if you are under that then just uploading the entire metadata is fine. 

In the app for uploading sample information, make sure to click the box "Ignore Warnings." Otherwise it says the cell runs successfully but it won't actually create the SampleSet if you have custom columns that trigger warnings.

Now you will have to link the workspace objects to samples, where you need to create another spreadsheet for linking the genome files to the samples in the metadata. For this, create a new narrative in KBase (such as sample-linkage) and create a code-cell. To create a code-cell press at the bottom right what looks like a command line. Insert this function, where `######` is your narrative ID which is in the URL link. 

```
# Create template for SampleSet Linking - Shareable
from biokbase.installed_clients.WorkspaceClient import Workspace
ws = Workspace("https://narrative.kbase.us/services/ws")
def print_sample_link_csv(wsid,objtype=None,typecol=True):
    '''
    Inputs:
        wsid - the workspace ID of the narrative you want a sample link template for
        objtype - Defaults "none" which returns all objects including the Narrative object. 
            Include an object type to select only that type.
        typecol - prints the type of the object to the csv, if you want to include all types so you can filter it later. 
            If you include this column, delete it from the final CSV before uploading to link samples.
    Outputs:
        Prints the text of a formatted CSV with all the information for all the objects you have selected except the sample name.
    Usage:
        Run this function for the Narrative that you want sample links for, then copy the output text to a text editor outside of KBase. 
        Fill in the "sample_name" column however is easiest for you, then upload as a csv/tsv/excel file and run "Batch link Workspace Objects to Samples" to perform the linking. 
        '''
    objects = ws.list_objects(params={"ids":[int(wsid)],"type":objtype})
    if typecol==True:
        csvtext = 'object_name,sample_name,obj. ref,DELETE_COLUMN-object_type'
    else:
        csvtext = 'object_name,sample_name,obj. ref'
    print(csvtext)
    for obj in objects:
        wsid = obj[6]
        objnum = obj[0]
        ver = obj[4]
        ref = f"{wsid}/{objnum}/{ver}"
        name = obj[1]
        if typecol==True:
            objtype = obj[2]
            line = f"{name},{name},{ref},{objtype}"
        else:
            line = f"{name},{name},{ref}"
        print(line)
        csvtext+=line
    return csvtext

csv = print_sample_link_csv(######)
```

The code-cell will output something that looks like this: 
```
object_name,sample_name,obj. ref,DELETE_COLUMN-object_type
Narrative.1748638437022,Narrative.1748638437022,218406/1/20,KBaseNarrative.Narrative-4.0
GCF_001434985.1_ASM143498v1_genomic,GCF_001434985.1_ASM143498v1_genomic,218406/2/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_003946605.1_ASM394660v1_genomic,GCF_003946605.1_ASM394660v1_genomic,218406/3/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_001436135.1_ASM143613v1_genomic,GCF_001436135.1_ASM143613v1_genomic,218406/4/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_003946555.1_ASM394655v1_genomic,GCF_003946555.1_ASM394655v1_genomic,218406/5/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_900618535.1_TH39_genomic,GCF_900618535.1_TH39_genomic,218406/6/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_010306265.1_ASM1030626v1_genomic,GCF_010306265.1_ASM1030626v1_genomic,218406/7/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_011067265.1_ASM1106726v1_genomic,GCF_011067265.1_ASM1106726v1_genomic,218406/8/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_001437055.1_ASM143705v1_genomic,GCF_001437055.1_ASM143705v1_genomic,218406/9/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_023091965.1_ASM2309196v1_genomic,GCF_023091965.1_ASM2309196v1_genomic,218406/10/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_003573535.1_ASM357353v1_genomic,GCF_003573535.1_ASM357353v1_genomic,218406/11/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_006007945.1_ASM600794v1_genomic,GCF_006007945.1_ASM600794v1_genomic,218406/22/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_003946305.1_ASM394630v1_genomic,GCF_003946305.1_ASM394630v1_genomic,218406/23/1,KBaseGenomeAnnotations.Assembly-5.1
GCA_003325435.1_Razy_CA_genomic,GCA_003325435.1_Razy_CA_genomic,218406/24/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_003946405.1_ASM394640v1_genomic,GCF_003946405.1_ASM394640v1_genomic,218406/25/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_001702115.1_ASM170211v1_genomic,GCF_001702115.1_ASM170211v1_genomic,218406/26/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_005405165.1_ASM540516v1_genomic,GCF_005405165.1_ASM540516v1_genomic,218406/27/1,KBaseGenomeAnnotations.Assembly-5.1
GCF_000420365.1_ASM42036v1_genomic,GCF_000420365.1_ASM42036v1_genomic,218406/28/1,KBaseGenomeAnnotations.Assembly-5.1
```

Note that for our specific purposes the filename of the assembly without the `.fa` extension is the name/sample name in the metadata, so the code-cell has been written to output the same string for `object_name` and `sample_name`. This might not be the case for your uploads. Save this code output as a CSV file. Upload it to the workspace. 

Now launch the Batch Link Workspace Objects to Samples app. 

## Add Documentation to the Narrative

## Create a Narrative with Subsets of Genomes