# PyMOL Interface Visualization Script for c7riboseqorf18_partner_01_SERPINE2.pdb
load C:\Users\Lucas.Escosteguy\documentos\microproteinproject\structures\pdbs\c7riboseqorf18_partner_01_SERPINE2.pdb, complex
hide everything, complex
show cartoon, complex
color marine, chain B
color forest, chain A

# Select and show interface residues
select if_uprot, chain A and resi 7+9+10+11+12+13+14+15+18+22+24+25+27+28+29+30+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52
select if_partner, chain B and resi 54+57+58+61+65+157+161+164+165+168+169+173+175+178+179+180+181+182+183+184+185+186+187+188+189+190+191+192+193+194+196+197+198+199+201+220+223+228+229+230+240+246+253+257+258+259+260+261+262+287+288+289+290+304+319+324+327+332+333+334+335+336+337+338+339+340+341+342+343+344+345+346+347+348+349+350+351+352+354+358+359+360+361+362+363+364+365+366+367+368+369+370+371+372+373+374+375+379+381+391

show sticks, if_uprot
show sticks, if_partner
color tv_yellow, if_uprot
color tv_orange, if_partner

# Highlight polar contacts / H-bonds
distance hbonds, chain A, chain B, 3.5, mode=2
color cyan, hbonds

set cartoon_transparency, 0.2
zoom if_uprot or if_partner
