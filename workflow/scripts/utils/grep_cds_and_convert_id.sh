for f in GFFs/* 
	do
		echo $f
		f2=${f##*/}
		echo $f2
		grep -w 'CDS' ${f} > GFFs_CDS/${f2}_CDS
	done

for f in GFFs_CDS/*
	do
		echo $f
		f2=${f##*/}
		echo $f2
		./gff_vbid_2_odbid.py ${f} GFFs_CDS_ODB/${f2}_ODB
	done

