f2py -c --noopt --debug --f77flags="-std=legacy -fbounds-check -fdefault-real-8" initialization.pyf fcclat.f bcclat.f nncal.f fccbnd.f bccbfe.f mmcal.f equiv.f
f2py -c --noopt --debug --f77flags="-std=legacy -fbounds-check -fdefault-real-8" slaterkoster.pyf slkode.f skdd.f selfd.f
f2py -c --noopt --debug --debug-capi --f77flags="-std=legacy -fbounds-check -fdefault-real-8" recal.pyf recal.f hop.f
f2py -c --noopt --debug --debug-capi --f77flags="-std=legacy -fbounds-check -fdefault-real-8" coeffproc.pyf recsum.f recqd.f recrts.f recwt.f cfgpgn.f cfgen.f denqd.f dencrs.f plyval.f termgn.f fenval.f bndref.f bndcrd.f bndest.f bndwt.f wtmin.f numd.f numc.f taban.f
