# homology
test_homology.ipynbがメインのソースコードです。
三角形分割による単体的複体Kをlist形式でcalc_HomologyGroupListに入力すると、
n次元ホモロジー群の基底をKのn単体を与えられた順に基底としたときの係数として出力し、さらに各基底に対するねじれ率torも取得されます。
計算の概要とtest_homology.ipynbにある例はhomology.pdfに記載しております。

計算の際に使用したスミス標準化のコードはhttps://qiita.com/yuji0001/items/64dc97cd4dcebf83d0a8 で紹介されていたものを使用しております。    

The main code is test_homology.ipynb.  
Given simplecial complex K obtained by triangulation, you get the pair of the basis set & the torsion of each basis vector for each n-th homology group.
The input K is given by the 2D list of the nodes corresponding with the triangulation.   
The output of the basis is the coefficiants of the n-th simplexies in the same order in K.
