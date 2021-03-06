import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


A = 'MTELKAKGPRAPHVAGGPPSP-EVGSPLLCRPAAGPFPGSQTSDTLPEVSAIPISLDGLLFPRPCQGQDPSDEKTQDQQSLSDVEGAYSRAEATRGAGGSSSSPPEKDSGLLDSVLDTLLAPSGPGQSQPSPPACEVTSSWCLFGPELPEDPPAAPATQRVLSPLMSRSGCKVGDSSGTAAAHKVLPRGLSPARQLLLPASESPHWSGAPVKPSPQAAAVEVEEEDGSESEESAGPLLKGKPRALGGAAAGGGAAAVPPGAAAGGVALVPKEDSRFSAPRVALVEQDAPMAPGRSPLATTVMDFIHVPILPLNHALLAARTRQLLEDESYDGGAGAASAFAPPRSSPCASSTPVAVGDFPDCAYPPDAEPKDDAYPLYSDFQPPALKIKEEEEGAEASARSPRSYLVAGANPAAFPDFPLG--PPPPLPPRATPSRPGEAAVTAAPASASVSSASSSGSTLECILYKAEGAPPQQGPFAPPPCKAPGASGCLLPRDGLPSTSASA-AAAGAAPALYPALGLNGLPQLGYQAAVLKEGLPQVYPPYLNYLRPDSEASQSPQYSFESLPQKICLICGDEASGCHYGVLTCGSCKVFFKRAMEGQHNYLCAGRNDCIVDKIRRKNCPACRLRKCCQAGMVLGGRKFKKFNKVRVVRALDAVALPQPVGVPNESQALSQRFTFSPGQDIQLIPPLINLLMSIEPDVIYAGHDNTKPDTSSSLLTSLNQLGERQLLSVVKWSKSLPGFRNLHIDDQITLIQYSWMSLMVFGLGWRSYKHVSGQMLYFAPDLILNEQRMKESSFYSLCLTMWQIPQEFVKLQVSQEEFLCMKVLLLLNTIPLEGLRSQTQFEEMRSSYIRELIKAIGLRQKGVVSSSQRFYQLTKLLDNLHDLVKQLHLYCLNTFIQSRALSVEFPEMMSEVIAAQLPKILAGMVKPLLFHKK'
arrayA = np.genfromtxt('P01624_sp_P06401_PRGR_HUMAN.out')
listA = []
i = 1
for letter in A:
    if letter == "-":
        listA.append(np.nan)
    else:
        listA.append(float(arrayA[i, 3]))
        i = i+1

B = 'MTELKSKGPRAPHVAGGPPSP-EVGSPLLCRPAAGPFQGSQTSDTLPEVSAIPISLDGLLFPRLCQGQDPPDKKTQNQQSLSDVEGAYSRAEATRGTGGSSSRPPEKDSGLLDSVLDTLLAPSGPGQSQPSPPACEVTSSWCLFGPELPEDPPAAPATQRVLSPLMSRSGGKTEDSSGTAAAHKVLPRGLSPSRQLLLPTSGSPHWSGAPVKPSPQPTAVEVEEEDGSESEDSAGPLLKGKSRVLGGAAAGGGAAAVPPGAAAGGVGLVPKEDSRFSAPRVALVEQDAPMAPGRSPLATTMMDFIHVPIVPLNHALLAARTRQLLEDESYDGGAGAASAFAPPQSSPSASSTPVAVGDFPDCAYPPDAEPKDNAYPLYGDFQPLALKIKEEEEGAEASARSPGSYLVAGANPAAFPDFPLG--PPPQLPPRAPPSRPGEAAVTAAPASASVSSASSPGSTLECILYKAEGALPQQGQFAPPPCKAPGAGGCLLPRDGLPSTSASAAAAAGAAPTLYPALGLNGLPQLGYQAAVLKEGLQQVYPPYLNYLRPDSEASQSPQYSFESLPQKICLICGDEASGCHYGVLTCGSCKVFFKRAMEGQHNYLCAGRNDCIVDKIRRKNCPACRLRKCCQAGMVLGGRKFKKFNKVRVMRALDAVALPQPVGIPNESQVLSQRFTFSPGQDIQLIPPLIKLLMSIEPDVIYAGHDNSKPDTSSSLLTSLNQLGERQLLSVVKWSKSLPGFRNLHIDDQITLIQYSWMSLMVFGLGWRSYKHVSGQMLYFAPDLILNEQRMKESSFYSLCLTMWQIPQEFVKLQVSQEEFLCMKVLLLLNTIPLEGLRSQTQFEEMRSSYIRELIKAIGLRQKGVVSSSQRFYQLTKLLDNLHDLVKQLHLYCLNTFIQSRALSVEFPEMMSEVIAAQLPKILAGMVKPLLFHKK'
arrayB = np.genfromtxt('P01624_sp_A7X8D2_PRGR_COLGU.out')
listB = []
i = 1
for letter in B:
    if letter == "-":
        listB.append(np.nan)
    else:
        listB.append(float(arrayB[i, 3]))
        i = i+1

C = 'MTELKAKGPRAPHVAGGPPSP-EVGSPLLCRPAAGPFPGSQTSDTLPEVSAIPISLDGLLFPRPCQGQDPSNEKTQDQQSLSDVEGAYSRAEATRGAGGSSSSPPEKDSGLLDSVLDTLLAPSGPGQSQPSPPACEVTSSWCLFGPELPEDPPAAPATQGVLSPLMSRSGCKAGDSSGTAAAHKVLPRGLSPSRQLLLPASGSPHWSGAPVKPSPQPAAVEVEEEDGSESEESAGPLLKGKPRALGGAAAGGGAAAVPPGAAAGGVALVPKEDSRFSAPRVALVEQDAPMAPGRSPLATTMMDFIHVPILPLNHALLAARTRQLLEDESYDGGAGAASAFAPPRSSPSASSTPVAVGDFPDCAYPPDAEPKDDAYPLYSDFQPPALKIKEEEEGAEASARSPRSYLVAGANPAAFPDFPLG--PPPPLPPRAPPSRPGEAAVTAAPASASVSSASSSGSTLECILYKAEGAPPQQGPFAPPPCKAPGASGCLLPRDGLPSTSASA-AAAGAAPALYPALGLSGLPQLGYQATVLKEGLPQVYPPYLNYLRPDSEASQSPQYSFESLPQKICLICGDEASGCHYGVLTCGSCKVFFKRAMEGQHNYLCAGRNDCIVDKIRRKNCPACRLRKCCQAGMVLGGRKFKKFNKVRVVRALDAVALPQPVGIPNESQALSQRFSFSPGQDIQLIPPLINLLMSIEPDVIYAGHDNTKPDTSSSLLTSLNQLGERQLLSVVKWSKSLPGFRNLHIDDQITLIQYSWMSLMVFGLGWRSYKHVSGQMLYFAPDLILNEQRMKESSFYSLCLTMWQIPQEFVKLQVSQEEFLCMKVLLLLNTIPLEGLRSQTQFEEMRSSYIRELIKAIGLRQKGVVSSSQRFYQLTKLLDNLHDLVKQLHLYCLNTFIQSRALSVEFPEMMSEVIAAQLPKILAGMVKPLLFHKK'
arrayC = np.genfromtxt('alpha_P0624_sp_A7X8B5_PRGR_PANPA.out')
listC = []
i = 1
for letter in C:
    if letter == "-":
        listC.append(np.nan)
    else:
        listC.append(float(arrayC[i, 3]))
        i = i+1

D = 'MTELKAKGPRAPHVAGGPPSP-EVGSPLLCRPAAGQFPGSQTSDTLPEVSAIPISLDGLLFPRPCQGQDPSYEKTQDQQSLSDVEGAYSRAEATRGAGGSSSSPPEKESGLLDSVLDTLLAPSGPRQSQPSPPACEVTSSWSLFGPELPEDPPAAPATQGVLSPLMSRSGGKAGDSSGTAAAHKVLPQGLSPSRQLLLPASGSPHWSGAPVKPSPQPAAVEVEEEDGSESEDSAGPLLKGKPRALGGAAAG-GAAAVPPGAAAGGVALVPKEDSRFSAPRVALVEQDAPMAPGRSPLATTVMDFIHVPILPLNHALLAARTRQLLEDENYDGGAGAASAFAPPRSSPSASSTPVAVGDFPDCAYPPDVEPKDDAYPLYGDFQPPALKIKEEEEGAEASARTPRSYLVAGANPAAFPDFPLG--PPPPLPPRAPPSRPGEAAVTAAPASASVSSASSSGSTLECILYKAEGAPPQQGPFAPPPSKAPGAGGCLPPRDGLPSTAASA-SAAGAAPALYPALRLNGLPQLGYQAAVLKEGLPQVYPPYLNYLRPDSEASQSPQYSFESLPQKICLICGDEASGCHYGVLTCGSCKVFFKRAMEGQHNYLCAGRNDCIVDKIRRKNCPACRLRKCCQAGMVLGGRKFKKFNKVRVVRALDAVALPQPVGIPNESQVLSQRITFSPGQDIQLIPPLINLLMSIEPDVIYAGHDNTKPDTSSSLLTSLNQLGERQLLSVVKWSKSLPGFRNLHIDDQITLIQYSWMSLMVFGLGWRSYKHVSGQMLYFAPDLILNEQRMKESSFYSLCLTMWQIPQEFVKLQVSQEEFLCMKVLLLLNTIPLEGLRSQTQFEEMRASYIRELIKAIGLRQKGVVSSSQRFYQLTKLLDNLHDLVKQLHLYCLNTFIQSRALSVEFPEMMSEVIAAQLPKILAGMVKPLLFHKK'
arrayD = np.genfromtxt('alpha_P0624_sp_A7X8C2_PRGR_HYLLA.out')
listD = []
i = 1
for letter in D:
    if letter == "-":
        listD.append(np.nan)
    else:
        listD.append(float(arrayD[i, 3]))
        i = i+1

df = pd.DataFrame(data={'P06401': listA,
                        'A7X8D2': listB,
                        'A7X8B7': listC,
                        'A7X8C7': listD})

mask = df.isnull()
ax = sns.heatmap(df, mask=mask)
ax.invert_yaxis()
plt.show()