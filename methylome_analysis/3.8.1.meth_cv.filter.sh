#awk '$3 <= 16573692 && $4 > 1 && $5 > 100 && $5 < 1000' CG.CV.bedGraph > CG.CV.filter.bedGraph
#awk '$3 <= 16573692 && $4 > 0.8 && $5 > 100 && $5 < 1000' CHG.CV.bedGraph > CHG.CV.filter.bedGraph
#awk '$3 <= 16573692 && $4 > 1 && $5 > 100 && $5 < 1000' CHH.CV.bedGraph > CHH.CV.filter.bedGraph
#awk '$3 <= 16573692 && $5 > 100 && $5 < 1000' CG.CV.bedGraph > CG.CV.filter.bedGraph
#awk '$3 <= 16573692 && $5 > 100 && $5 < 1000' CHG.CV.bedGraph > CHG.CV.filter.bedGraph
#awk '$3 <= 16573692 && $5 > 100 && $5 < 1000' CHH.CV.bedGraph > CHH.CV.filter.bedGraph
awk '$4 > 0' CG.CV.bedGraph > CG.CV.before.bedGraph
awk '$4 > 0' CHG.CV.bedGraph > CHG.CV.before.bedGraph
awk '$4 > 0' CHH.CV.bedGraph > CHH.CV.before.bedGraph

# awk '$5 > 100' CG.CV.bedGraph > CG.CV.filter.bedGraph
# awk '$5 > 100' CHG.CV.bedGraph > CHG.CV.filter.bedGraph
# awk '$5 > 100' CHH.CV.bedGraph > CHH.CV.filter.bedGraph
