# PyMOL Interface Visualization Script for cXriboseqorf59_partner_02_KAT5.pdb
load C:\Users\Lucas.Escosteguy\documentos\microproteinproject\structures\pdbs\cXriboseqorf59_partner_02_KAT5.pdb, complex
hide everything, complex
show cartoon, complex
color marine, chain B
color forest, chain A

# Select and show interface residues
select if_uprot, chain A and resi 6+9+10+12+13+14+16+17+18+19+20+21+22+23+24+25+26+27+28
select if_partner, chain B and resi 426+430+436+437+438+439+440+441+442+452+453+456+460+465+472+473+474+475+476+477+478+481

show sticks, if_uprot
show sticks, if_partner
color tv_yellow, if_uprot
color tv_orange, if_partner

# Highlight polar contacts / H-bonds
distance hbonds, chain A, chain B, 3.5, mode=2
color cyan, hbonds

set cartoon_transparency, 0.2
zoom if_uprot or if_partner
