# 2Link-Manipulator-Simulation
2リンクマニピュレータの動力学シミュレーション。
このシミュレーターでシンプルな制御システムを身につけることを目的としている。
そのため、外部のライブラリは使用していない。

![figure](https://cdn-ak.f.st-hatena.com/images/fotolife/s/sin6rai/20080927/20080927001629.jpg)
質力　　　 m1, m2
リンク長　　 L1, L2
リンク中心長 Lc1, Lc2　
イナーシャ 　 J1, J2
重力加速度 　 g
関節トルク τ1,τ2

𝕄(𝕢)𝕢¨+ℂ(𝕢,𝕢˙)+𝔾(𝕢)=τ

M11=J1+J2+m1Lc21+m2L21+m2Lc22+2m2L1Lc2cos(θ2)
M12=J2+m2Lc22+m2L1Lc2cos(θ2)
M21=M12
M22=J2+m2Lc22
C1=−m2L1Lc2sin(θ2)(2θ1˙θ2˙+θ2˙2)
C2=m2θ1˙2L1Lc2sin(θ2)
G1=−m1gLc1sin(θ1)−m2g(L1sin(θ1)+Lc2sin(θ1+θ2))
G2=−m2gLc2sin(θ1+θ2)


コントローラは目標角度に追従するようなPDコントローラにすると、
τ1=cp1(dv1−θ1)−cd1θ1˙
τ2=cp2(dv2−θ2)−cd2θ2˙
cp1,cp2,cd1,cd2 : コントローラ係数
dv1,dv2 : 目標角度

![figure2](https://cdn-ak.f.st-hatena.com/images/fotolife/s/sin6rai/20080930/20080930162420.jpg)
システム系入力  :  θd˙ 、 θd
アームモデル入力  :  コントローラ出力τ
システム系出力  :  θ˙ 、 θ
