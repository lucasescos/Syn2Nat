# PyMOL Interface Visualization Script for c10norep31_partner_05_PPP1CA.pdb
load C:\Users\Lucas.Escosteguy\documentos\microproteinproject\structures\pdbs\c10norep31_partner_05_PPP1CA.pdb, complex
hide everything, complex
show cartoon, complex
color marine, chain B
color forest, chain A

# Select and show interface residues
select if_uprot, chain A and resi 3+4+6+7+8+9+12+19+20+21+22+23+24+25+26+27+28+29
select if_partner, chain B and resi 48+49+50+53+54+55+56+78+116+119+166+167+168+169+242+243+253+257+261+279+283+288+289+290+291+292+293+294+295+296+297

show sticks, if_uprot
show sticks, if_partner
color tv_yellow, if_uprot
color tv_orange, if_partner

# Highlight polar contacts / H-bonds
distance hbonds, chain A, chain B, 3.5, mode=2
color cyan, hbonds

set cartoon_transparency, 0.2
zoom if_uprot or if_partner
