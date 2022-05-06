# Uniref2Lineage
A python script to support UniprotKB mapping from uniref accessions to uniprot accessions and uniprot metadata, such as lineages go terms and enzyme classification numbers.

<br>
<br>Input paramaters are:
<br>accs: a list of your accessions 
<br>input_data_type: the type of accession that you input (https://www.uniprot.org/help/api_idmapping)
<br>columns: data columns required in the final output (https://www.uniprot.org/help/uniprotkb_column_names)



<br>
<br>The script works by first mapping Uniref accessions to UniprotKB accessions.
<br>Then these accessions are mapped again to UniprotKB with specified colunms. 
<br>Since the script uses multithreading, please do not change too much the parameters "accession_batch" "batch_no", as increasing these could crash the Uniprot website (it has happened to me!)
