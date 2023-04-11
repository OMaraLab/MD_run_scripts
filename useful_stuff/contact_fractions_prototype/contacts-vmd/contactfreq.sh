atomselect macro GlyT2 {resname GLY ALA SER THR CYS VAL LEU ILE MET PRO PHE TYR TRP ASP GLU ASN GLN HIS LYS ARG}
atomselect macro ligand {resname BB4R BB4S C06R C06S E25R E25S G07R G07S H03R H03S}


contactFreq {GlyT2} {ligand} 4 0 contact-ligand.out
