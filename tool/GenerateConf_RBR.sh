#!/bin/bash 
set -o noclobber
for i in {5351..6372}
#for i in {5351,5352}
do
	file="./param/conf_RBR/analyzer_e42RBR_Run0$i.conf"
	echo "Generating $file"
	pK18=`grep pK18 /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$i | cut -c17-28`
	Kurama_F=`grep KURAMA_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$i | cut -c18-28`
	SHS_F=`grep SHS_F /group/had/sks/E42/JPARC2021May/K18MagnetParam/K18MagParam_0$i | cut -c18-28`
	echo "#">>$file
	echo "#  Generated Conf file for Run By Run PHC: run 0$i ">>$file
	echo "#">>$file
	echo "">>$file
	echo "UNPACK:		param/UNPACK/unpacker.xml">>$file
	echo "DIGIT:		param/DIGIT/digit_data_e42_20210423.xml	">>$file 
	echo "CMAP:		param/CMAP/channel_map_e42_20210423.xml">>$file 
	echo "DCGEO:		param/DCGEO/DCGeomParam_E42_20210606">>$file 
	echo "DCDRFT:		param/DCDRFT/DCDriftParam_05047">>$file 
	echo "DCTDC:		param/DCTDC/DCTdcCalib_05047">>$file 
	echo "HDPRM:		param/HDPRM/HodoParam_E42_20210708">>$file 
	echo "HDPHC:		param/HDPHC_RBR/HodoPHCParam_run0$i">>$file
	echo "USER:		param/USER/UserParam_E42_2021May_default">>$file 
	echo "#TPCPRM:		param/TPCPRM/TPCParam_0">>$file 
	echo "TPCPRM:		param/TPCPRM/TPCParam_Defocus2_ON">>$file 
	echo "TPCPOS:		param/TPCPOS/TPCPositionCorrectionMap_0">>$file 
	echo "#TPCPOS:		param/TPCPOS/TPCPositionCorrectionMapX_defocus_paramy_DefocusON_ym60p60_mesh15">>$file 
	echo "BH2FLT:		param/BH2FLT/BH2Filter.param">>$file 
	echo "BH1MTH:		param/BH1MTH/BH1Matching.param">>$file 
	echo "">>$file
	echo "FLDCALC:	0.749779">>$file 
	echo "FLDNMR:		$Kurama_F">>$file 
	echo "FLDMAP:		param/KURAMA/KuramaFieldMap_20210526">>$file 
	echo "">>$file
	echo "HSFLDCALC:	0.90607">>$file
	echo "HSFLDHALL:	$SHS_F">>$file
	echo "HSFLDMAP:	param/KURAMA/ShsFieldMap_20210526">>$file
	echo "">>$file
	echo "#PK18:		0.5970048">>$file 
	echo "#K18TM:		param/K18TM/K18MatrixParam_0.6_p_05311">>$file 
	echo "PK18:		$pK18">>$file 
	echo "K18TM:		param/K18TM/K18MatrixParam_2019Feb_plus">>$file 
done
