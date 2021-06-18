## Script to convert a subsitutions per site tree
## into an ultra-metric time tree.

## Loading required library
library(ape)

## Setting tree directory path
setwd("~/evol-feat-snakemake/workflow/input/6656/")

## We read in the PHUMU+RPROC rooted tree, 'RO'.
## The tree was rerooted with Newick utilities, terminal command:
## nw_reroot Insecta_moz2_tree.pep PHUMU RPROC | nw_order -> Insecta_moz2_tree_RO.pep

## Reading tree file
molecular_tree <- read.tree("arthropoda_rooted.tree")

## Or reading characters from file directly in script
#molecular_tree <- read.tree(text="((PHUMU:0.52318811260609177083,RPROC:0.57574984435598719479)100:0.0185631,((AMELL:0.13420825900633098771,LHUMI:0.14693184769600894923)100:0.31233543682593423174,(((((((((((DERE1:0.01021522022447337620,DYAK1:0.01106527559303413337)100:0.00401176445984145497,((DSEC1:0.00814300186400131361,DSIM1:0.00962408213703389546)100:0.00552762487274174007,DMEL5:0.00749867382839301908)100:0.00919038874702339005)100:0.03760133861117224280,DANA1:0.04974961170084617190)100:0.02382661163059059811,(DPER1:0.01309792385561942343,DPSE3:0.00224195353017842946)100:0.06136856290688777471)100:0.02482842762886991461,DWIL1:0.08907899707058640482)100:0.01880645625904252966,((DMOJ1:0.04907090912886290929,DVIR1:0.03106631254492941280)100:0.01423643328465117562,DGRI1:0.05642105389690358491)100:0.04027582044845061382)100:0.14161192708217770764,GMORY:0.21316868653477499818)100:0.25289136828831959569,(LLONJ:0.08531352032017081033,PPAPI:0.10223252622284993707)100:0.24871346684319486919)0:0.03670333192798320293,(((AALBS:0.01705456225914456572,ADARC:0.01773903381715946984)100:0.07396308046482585885,((AATRE:0.03438791310038371435,ASINS:0.03420038554527971703)100:0.03643876535770280978,((ADIRW:0.02197714162279932118,AFARF:0.02475040438528573750)100:0.02388059772103924908,((((ACULA:0.01698004526448793419,AMINM:0.01454013585743476646)100:0.00745115218901763830,AFUNF:0.01890321968740657846)100:0.00993182975318094956,((ASTEI:0.00152087762278565932,ASTES:0.00169765165208231599)100:0.01439242356033401497,AMACM:0.01922076860097963330)100:0.01333152049049606339)100:0.01215478545591066988,(((((AMELC:0.00616219782682916710,AMERM:0.00339193633835763869)88:0.00041900124007265100,(AARAD:0.00338699885452365627,AQUAS:0.00303440811615654244)98:0.00082941788173386294)90:0.00083922864181178430,AGAMP:0.00355552207143021350)100:0.02087080000156923271,ACHRA:0.02508470663698924435)100:0.00552075492727747107,AEPIE:0.02823599498816542486)100:0.01186417572110551322)100:0.01257086921689764061)100:0.03151476363212395160)100:0.01976105072198313176)100:0.09178695837509466549,(AAEGL:0.07695631746339641477,CPIPJ:0.09619389850005723164)100:0.06627564187558378228)100:0.16763633397699556626)100:0.23118713943721178050,(BMORI:0.16953288605733374261,DPLEX:0.16180401498822019613)100:0.48039750964988270354)100:0.07504104418349075156,TCAST:0.40678949602550951159)100:0.06822135232405733551)100:0.0556892);")

## Plotting tree
plot(molecular_tree)

## Checking if tree is ultra-metric
is.ultrametric(molecular_tree)

## Calibrating time from estimation of divergence between P. humanus
## and all other species (i.e. 366 million years ago)

nodes <- c(
  getMRCA(molecular_tree, tip = c("Strigamia_maritima","Drosophila_melanogaster") ),
  getMRCA(molecular_tree, tip = c("Apis_mellifera","Limulus_polyphemus") ),
  getMRCA(molecular_tree, tip = c("Apis_mellifera","Tribolium_castaneum") ),
  getMRCA(molecular_tree, tip = c("Apis_mellifera","Bombus_impatiens") ),
  getMRCA(molecular_tree, tip = c("Apis_mellifera","Megachile_rotundata") ),
  getMRCA(molecular_tree, tip = c("Glossina_morsitans","Anopheles_gambiae") ),
  getMRCA(molecular_tree, tip = c("Aedes_aegypti","Anopheles_gambiae") ),
  getMRCA(molecular_tree, tip = c("Calopteryx_splendens","Blattella_germanica") ),
  getMRCA(molecular_tree, tip = c("Anopheles_atroparvus","Anopheles_gambiae") ),
  getMRCA(molecular_tree, tip = c("Parasteatoda_tepidariorum","Strigamia_maritima") ),
  getMRCA(molecular_tree, tip = c("Danaus_plexippus","Pediculus_humanus") )
) 

calib <- makeChronosCalib(molecular_tree, node=nodes, age.min=c(583, 601, 325, 102, 110, 272, 190, 413, 52, 601, 358), age.max=c(583, 601, 325, 102, 110, 272, 190, 413, 52, 601, 358)+1)

## Converting substitutions per site branch lenght to time branch lengths
ultrametric_tree <- chronos(molecular_tree, calibration = calib, lambda = 1,
                            model = "discrete")

## Plotting tree
plot(ultrametric_tree)

## Checking if tree is ultra-metric
is.ultrametric(ultrametric_tree)

## ATTENTION: in some cases, the ultrametrization of the tree might round the 
## branchlengths to an extent in which " is.ultrametric() " will no longer
## recognize an ultrametric tree, due to the small rounding errors. If this happens,
## check again by setting a tolerance parameter (adjustable, it has to be much
## smaller than the differences between close branchlengths) as follows:
is.ultrametric(ultrametric_tree, tol = 0.0001)


## Writing ultrametric tree into a new file
write.tree(ultrametric_tree, file = "arthropoda_rooted_time.tree")


## Additionally writing a second tree with rounded branch lengths
ultrametric_tree_rounded <- roundBranches(ultrametric_tree, digits = 3)
is.ultrametric(ultrametric_tree_rounded)
write.tree(ultrametric_tree_rounded, file = "trees/Insecta_ultrametric_tree_rounded.pep")
