# **KRE_DOP: Key Residues Explored in Allostery of *δ* Opioid Receptor Using Elastic Network Model and Complex Network Model Combined with Perturbation Method.**
Authors are Lei Chen, Weikang Gong, Zhongjie Han, Wenxue Zhou, Shuang Yang, Chunhua Li

The performance process includes three steps: data preparation, data extraction, fcfGNMMD constructed and key residues explored

KRE_DOP uses the following dependencies:
- Matlab 2019b
- R 4.0.3 
- VMD

## The process of calculation:
### Data preparation
Molecular dynamics (MD) simulation trajectory was available in GPCRmd database, including “10714_trj_73.xtc” and “pdb0.pdb”.

File “10714_trj_73.xtc” has a large mempry, so you need to download it from the GPCRmd database.

The crystal structure of δ Opioid Receptor (DOP) was obtain in RCSB database called“4n6h_A.pdb”. 

These above data are in the folder “DATA”. 
### Data extraction
#### Step1:
Run VMD with“pdb0.pdb” and “10714_trj_73.xtc” as input file, then you can obtain “1074trj-cas.dcd”.

Run “mddata_R.R” and you can obtain Root-mean-square deviation (RMSD) of the trajectory. The MD simulation trajectory began to stabilize in 100-500ns (499-2499 frame).
#### Step2:
Run VMD with“pdb0.pdb” and “10714_trj_73.xtc” as input file, then you can obtain “1074trj-cas100-500ns.dcd” and “499CA.pdb”.
#### Step3:
Run ”md_R.R” and obtain the files called “md_msf.xlsx, md_cof.xlsx and md_cij.csv", respectively.

These above data and code are in the folder “MD_R” .
### fcfGNMMD constructed and key residues explored
Run “fcfGNMMD.m”, “DPR.m” and “compnetwork.m”.
The output is the value of PCC, the dissipated work, the dynamic correlations upon a site perturbed, degree and Z-score.

These above data and code are in the folder “KER-MATLAB”.

## Help 
For any questions, please contact us by chunhuali@bjut.edu.cn.
