# Examples

This page lists some examples on how to use the genie package.

## KRAS G12C

<div>
  <table>
    <tbody>
      <tr>
        <td><b>Program:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/kras_g12c.py">kras_k12c.py</a</td>
      </tr>
      <tr>
        <td><b>Result:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/KRAS_G12C_GENIE_15.0.xlsx">KRAS_G12C_GENIE_15.0.xlsx</a></td>
      </tr>
    </tbody>
  </table>
</div>

This program creates an Excel file with the frequencies of KRAS G12C mutations in non-small cell lung cancer (NSCLC) and in colorectal cancer (CRC).

Frequencies are provided for all NSCLC and for all CRC samples, as well as for some cancer subtypes:

* NSCLC: LUAD (lung adenocarcinoma), LUSC (lung squamous cell carcinoma), other, all
* CRC: COAD (colon adenocarcinoma), READ (rectal carcinoma), COADREAD (colorectal carcinoma), other, all

Frequencies are provided for each race (White, Black, Asian, other, all).

## HER2 (ERBB2) - top cancers by race

<div>
  <table>
    <tbody>
      <tr>
        <td><b>Program:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/her2.py">her2.py</a</td>
      </tr>
      <tr>
        <td><b>Result:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/HER2_GENIE_15.0_top_cancers.xlsx">HER2_GENIE_15.0_top_cancers.xlsx</a></td>
      </tr>
    </tbody>
  </table>
</div>

This program creates an Excel file with the frequencies of all amino acid variants found in the five cancer types with the most frequent HER2 mutations:

* bladder cancer
* small bowel cancer
* cervical cancer
* non-small cell lung cancer
* esophagogastric cancer

Separate frequencies are provided for the mosts frequent subtypes of these cancers. Furthermore, separate allele frequencies are provided for each race (White, Black, Asian, other, all).

## HER2 (ERBB2) - top cancers by TMB

<div>
  <table>
    <tbody>
      <tr>
        <td><b>Program:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/her2_with_tmb.py">her2_with_tmb.py</a</td>
      </tr>
      <tr>
        <td><b>Result:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/HER2_GENIE_15.0_top_cancers_with_tmb.xlsx">HER2_GENIE_15.0_top_cancers_with_tmb.xlsx</a></td>
      </tr>
    </tbody>
  </table>
</div>

This program creates an Excel file with the frequencies of all HER2 amino acid variants found in the five cancer types with the most frequent HER2 mutations, same as described in the previous section, but here separate frequencies are provided for each tumor mutation burdon class (_low_, _intermediate_, _high_) instead of for each race.

## HER2 (ERBB2) - all cancers

<div>
  <table>
    <tbody>
      <tr>
        <td><b>Program:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/her2_andi.py">her2_andi.py</a</td>
      </tr>
      <tr>
        <td><b>Result:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/HER2_GENIE_15.0_all_cancers.xlsx">HER2_GENIE_15.0_all_cancers.xlsx</a></td>
      </tr>
    </tbody>
  </table>
</div>

This program creates an Excel file with the frequencies of all HER2 amino acid variants that are not functionally deleterious  (such as frame shifts, stop gained or lost, ...) for each cancer type. In addition to amino acid level mutation frequencies, gene level frequencies are provided.

## HER2 (ERBB2) - clinical trial indications, different mutation sets

<div>
  <table>
    <tbody>
      <tr>
        <td><b>Program:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/her2_andi_v2.py">her2_andi_v2.py</a</td>
      </tr>
      <tr>
        <td><b>Result:</b></td>
        <td><a href="https://github.com/Bayer-Group/genie.parser/blob/main/examples/HER2_GENIE_15.0_selected_cancers_gene_level.xlsx">HER2_GENIE_15.0_selected_cancers_gene_level.xlsx</a></td>
      </tr>
    </tbody>
  </table>
</div>

This program creates an Excel file with the gene level mutation frequencies of HER2. Several definitions for the gene level were used. Each of these definitions is based on a different set of amino acid changes, and HER2 is classified as mutated if a sample has at least one of the mutations from the respective set. The different mutation sets were defined by the HER2 project team:

* **curated:**   variants that have been curated as activating HER2 mutations
* **phase3:**    variants that have been selected as eligible for inclusion in the NSCLC Phase 3 clinical study
* **tkd:**       variants that are located in the Tyrosine Kinase Domain of HER2; TKD is defined as exons 18 - 22; Pfam TKD model matches AAs 721 - 975
* **exon20ins:** variants that are an inframe exon20 insertion mutation
* **missense:**  variants that are a SNV missense mutation
