# %% markdown
# # rightMINI project data analysis log book
#
# ## Introduction
# The rightMINI analysis is completely based in the MINI project.
# The preprocessing was done in MINI. Refer to that documentation
# Steps taken here are:
#   1. Create the rhVOTC
#   2. Run batch analyze
#   3. Move data to the surface and go from there
#
# %%
# Primero hacer todos los imports y asegurar que chusca bien
import pandas as pd
import re
import os, shutil
import subprocess as sp
import glob as glob

# Function definitions
import glutil  # functions created by me in a different module.
join = os.path.join

# useful directories
print os.getcwd()
basedir = '/bcbl/home/public/Gari/MINI'
import platform
if platform.node()=='GariMBP': basedir = '/Users/gari/Documents/BCBL_PROJECTS/MINI'
DATAdir= join(basedir, 'DATA')
dicomdir=join(DATAdir, 'dicoms')
ANALYSISdir = join(basedir, 'ANALYSIS')
fMRI_SPMdir = join(ANALYSISdir, 'fMRI_SPM')
SUBJECTS_DIR = join(ANALYSISdir, 'freesurfer')# Freesurfer SUBJECTS_DIR
SUBJECTS_DIRacpc = join(ANALYSISdir, 'freesurferacpc')# Freesurfer acpc SUBJECTS_DIR
qMRIdir = join(ANALYSISdir, 'qMRI')
DWIdir = join(ANALYSISdir, 'DWI')
retdir = join(ANALYSISdir, 'ret')
os.chdir(basedir)
print os.getcwd()
# %%
# Leer todo el pedazo de archivo
SEGU = pd.read_csv('https://docs.google.com/spreadsheets/d/1LQnX8nsXr-lZQqOziO8qiiM_I1huGiCNxdUHCY7mybI/' +
                    'export?gid=0&format=csv')
# ANALYSISdir = "~/Desktop/ANALYSIS"
# SEGU = pd.read_excel(join(ANALYSISdir, 'SEGU.xlsx'))
GLMs = ['event', 'block']


print 'END OF CELL'
# %%
# SEGU.to_csv(join(ANALYSISdir,'tempSEGU3.csv'), sep='\t', encoding='utf-8')
# %%
SEGU.loc[0:49, ['MINID', 'artRep_event','artRepR1_block', 'artRepR2_block', 'artRepR1_triad', 'artRepR2_triad',
                'reconAll', 'reconAllOld', 'mrQnifti']]
# %%
SEGU.loc[50:99, ['MINID', 'artRep_event','artRepR1_block', 'artRepR2_block', 'artRepR1_triad', 'artRepR2_triad',
                'reconAll', 'reconAllOld', 'mrQnifti']]
#
# - PP1:Preprocessing 1: a:slice_timing + r:realignment + s:partial smoothing 4 + v:artifact repair + w:normalizacion (corregistracion a segmentacino  + normalizacion)
# - PP2:Preprocessing 2:QA:
#      - Crear carpeta en data que sea non_valid y dentro /motion y /otros /outlier en behav etc
#      - leer archivos artifact repair por RUN y hacer conteo de los valores que aparecen (dentro de cada run habra un archivo que se llama art_repaired con el num de slices, queremos el num). Sacar media de voluemnes corregidos por run y luego la media del sujeto a traves del run. decidir si se elimina o no: poner threshold comun: 15% 20% lo mas elegante 10%. Puede pasar que vovlamos a correr PP1 cambiando de 0.5 a 1mm por ejemplo.
#     - (Solo mirar los que artifact repair no haya detectado) ver que no tienen drift de mas de 3mm dentro de cada run. Eliminarlo si esta muy mal. Moverlo a /nonvalid/motion y marcarlo en el excel
# - PP3: s: smoothing final con 3mm, hay que cambiar en el experiment.m donde pone 4 pasar a 3 y el comand line se pasa solo el s.
# - PP4: QA
#      - mirar las imagenes check
#      - como estamos haciendo segmentacion 1 de cada 10 dara problemas con la normalizacon. Hay que revisar la mascara de la normalizacion:
#
#
#
# ## Group Analysis
# ### Preproc: arsvws
# - a: slice timing correction
# - r: realigment
# - s: 4mm smoothing
# - v: volume correction with artifact repair
# - w: normalization to the MNI space
# - s: 3mm smoothing, to make it 5 in total, since sqrt(4^2 + 3^2)= 5
#
#
#
#
# ## Individual Analysis
# ### Preproc: arsvct
# - a: slice timing correction
# - r: realigment
# - s: 4mm smoothing
# - v: volume correction with artifact repair
# - c: coregistration of the functional files to the anatomical
# - t: reslicing of the 2.5^3 data to 1^3 data, in the anatomical space. Now all info is in the anatomical acpc space.
#
#

# %% markdown

# APUNTAR AQUI TODOS LOS SUJETOS QUE NO TENIAN EL ESPERADO NUM DE ARCHIVOS Y LO QUE SE HIZO CON ELLOS
# (todos los fallos iniciales eran porque habia carpetas con run-s inacabados. Limpar esto para solo dejar aquellos que no tengan el mismo numero de dicoms y esto afecte a SPM y que quede aqui reflejado. )
#
# - S012 S012_MINIONCE_3940
#     - event
#         - ERROR: apuntar sujeto, no tiene 645 archivos como resto de event
#         - Run esta vacio
#         - temp esta vacio
#         - habia dos intentos EVENT y el primero estaba vacio, borrado y se vuelve a hacer mcconvert, pero solo para event
#         - `qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/event -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S012_MINIONCE_3940/EVENT_RUN1de1_ep2d_Mini_ipat2_TR2400_10gap_29"`
#
# - S013 S013_MINIONCE_3843
#     - triad
#         - ERROR: apuntar sujeto, no tiene 405 archivos como resto de triad
#         - `qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/triad -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S013_MINIONCE_3843/TRIAD_RUN1de2_ep2d_Mini_80"`
#     - block
#         - ERROR: apuntar sujeto, no tiene 575 archivos como resto de block
#         - Run1 estaba vacio, mismo problema que event arriba, varios intentos. Esto hay que limpiarlo antes de lanzar...
#         - `qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/block -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S013_MINIONCE_3843/BLOQUE_RUN1de2_ep2d_Mini_ipat2_TR2400_10gap_48"`
#
#
# NO VOY A HACERLOS TODOS A MANO: voy a eliminar los runs malos de la carpeta dicoms, y volvere a lanzar los sujetos enteros.
# Aqui abajo ire apuntando todos los que dieron problemas y se borraron runs a mano.
# PONER DO = 1 SOLO en los siguientes sujetos: 17 y 22
#
#
#
#
#
# - S017 S017_MINIONCE_4777 >>>>>>> S017_DAY1_MINITWICE_4777
#     - event
#         - ERROR: apuntar sujeto, no tiene 645 archivos como resto de event
#     - block
#         - ERROR: apuntar sujeto, no tiene 575 archivos como resto de block
#     - Se vuelve a hacer entero, dicoms arreglados
#     - Nota David: Día 1 EVENT ha quedado 1 estímulo fuera de la adquisición. En el bloque A entre bloque y retinotopía el SCANNER da error y reiniciamos el programa SYNGO. Hago las retinotopías sin localizer para ganar tiempo
#
# - S019 S019_MINIONCE_3740
#     - event
#         - ERROR: apuntar sujeto, no tiene 645 archivos como resto de event
#         - ARREGLAR A MANO: run1 29  ` qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/event -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S019_MINIONCE_3740/EVENT_RUN1de1_ep2d_Mini_ipat2_TR2400_10gap_29" `
#
# - S022 S022_MINIONCE_5276
#     - event
# ERROR: apuntar sujeto, no tiene 645 archivos como resto de event
#     - triad
#         - ERROR: apuntar sujeto, no tiene 405 archivos como resto de triad
#         - OJO: mirar el logs, abajo hemos visto que triad es problematico
#     - block
# ERROR: apuntar sujeto, no tiene 575 archivos como resto de block
#     - NOTA DAVID: EVENT 1 Error DirectX, TRIADS RUN 1 da error DirectX cuando quedaba 1 minuto. Reiniciamos PC y corremos sin errores. No podemos correr retinotopia, viene un segundo día (18/11)
#     - TRIAD: da error de nuevo, por lo qu edice David me imagino, muevo archivos a mano pero ojo con esto al calcular
#
# Los siguientes a mano tb
# - S023 S023_MINIONCE_5280
#     - block
#     - ERROR: apuntar sujeto, no tiene 575 archivos como resto de block
#     - A mano, solo tenia dos bloques repetidos: ` qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/block -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S023_MINIONCE_5280/BLOQUE_RUN1de2_ep2d_Mini_ipat2_TR2400_10gap_29"`
#
# - S024 S024_MINIONCE_5286
#     - block
#         - ERROR: apuntar sujeto, no tiene 575 archivos como resto de block
#         - `  qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/block -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S024_MINIONCE_5286/BLOQUE_RUN1de2_ep2d_Mini_ipat2_TR2400_10gap_28" `
#
# - S025 S025_MINIONCE_3096
#     - event
#         - ERROR: apuntar sujeto, no tiene 645 archivos como resto de event
#         - ` qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/event -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S025_MINIONCE_3096/EVENT_RUN1de1_ep2d_Mini_ipat2_TR2400_10gap_8"  `
#
# - S027 S027_MINIONCE_2367
#     - block
#         - ERROR: apuntar sujeto, no tiene 575 archivos como resto de block
#         - run1 7:   `  qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/block -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S027_MINIONCE_2367/BLOQUE_RUN1de2_ep2d_Mini_ipat2_TR2400_10gap_7"  `
#
# - S029 S029_MINIONCE_1459
#     - event
#         - ERROR: apuntar sujeto, no tiene 645 archivos como resto de event
#         - Aqui el problema era que habia repetidos pq vino el lehendakari e hicieron el paripe, habia un log que se llamaba fake2 tb, ojo por si vuelve a aparecer, creo que ya esta borrado:
#         - run1 27  ` qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/event -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S029_MINIONCE_1459/EVENT_RUN1de1_ep2d_Mini_ipat2_TR2400_10gap_27"  `  NO LO LANZO, ya esta hecho, revisar que pasa
#         - Faltaba un dicom, se vuelve a exportar de osirix y volver a convertir.
#
# - S030 S030_MINIONCE_1616
#     - block
#         - ERROR: apuntar sujeto, no tiene 575 archivos como resto de block
#         - run1 55   `  qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/block -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S030_MINIONCE_1616/BLOQUE_RUN1de2_ep2d_Mini_ipat2_TR2400_10gap_55"  `
#
# - S031 S031_MINIONCE_1144
#     - event
#         - ERROR: apuntar sujeto, no tiene 645 archivos como resto de event
#         - run1 39  ` qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/event -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S031_MINIONCE_1144/EVENT_RUN1de1_ep2d_Mini_ipat2_TR2400_10gap_39" `
#
# - S032 S032_MINIONCE_5431
#     - event
#         - ERROR: apuntar sujeto, no tiene 645 archivos como resto de event
#         - run1 31  ` qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/event -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S032_MINIONCE_5431/EVENT_RUN1de1_ep2d_Mini_ipat2_TR2400_10gap_31"  `
#
#
# - S048 S008_DAY1_MINITWICE_2511
# event
# triad
# block
# ERROR: apuntar, 11 f_s, no tiene 575 archivos como resto de block
#
#     - LAnamos solo el bloc, haia unos pocos dicoms dentro que ha hecho que se lie el asunto:
#     - run
#     - REVISAR ESTE SUJETO A MANO
# qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/block -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S008_DAY1_MINITWICE_2511/BLOQUE_RUN1de2_ep2d_Mini_ipat2_TR2400_10gap"
#
#     - MOVE: src: IM-0007-0001.dcm dest: ../block2normalizados/IM-0007-0001.dcm
# qsub -v DIR=/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM $mySH/RunPython.sh "mcverter -o /bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/temp_func_files/block -f nifti -q -j -x -F f_ /bcbl/home/public/Gari/MINI/DATA/dicoms/S008_DAY1_MINITWICE_2511/BLOQUE_RUN2de2_ep2d_Mini_ipat2_TR2400_10gap"
#
#     - Arreglado
#
#
# - S078 S008_DAY2_MINITWICE_2511
# event
# triad
# block
# ERROR: apuntar, 413 f_s, no tiene 575 archivos como resto de block
#     - Arreglado, lo he hecho a mano

# %% markdown
# ## LOG files
#
# Hacer los logs usando el script que hice hace tiempo: sReadMRIlogs.py
#
# ```bash
# python $c/MINI/sReadMRIlogs.py -i /bcbl/home/public/Gari/MINI/DATA/logs/MRI/block -o FORSPM -r 2 -t Block -p mrilog
# python $c/MINI/sReadMRIlogs.py -i /bcbl/home/public/Gari/MINI/DATA/logs/MRI/event -o FORSPM -r 1 -t EVENT -p mrilog
# python $c/MINI/sReadMRIlogs.py -i /bcbl/home/public/Gari/MINI/DATA/logs/MRI/triad -o FORSPM -r 2 -t Triad -p mrilog
# ```
#
# - Estas las hemos hecho directmaente en el command line.
# - Ahora mover todos los archivos convertidos a converted y los csv de FORSPM a Analysis.
# Validar que tenga la longitud adecuada y apuntar aqui lo que tengan problemas como referencia:
#
# **BLOCK**
# - No esta creando el S007_DAY2_MINITWICE, ver que pasa:
#     - El problema era que tenian un 1 al final del nombre del archivo
#
# **EVENT**
# - S004_DAY1_MINITWICE_3011: se hizo dos veces y las dos veces fallo en el mismo sitio, aunque llego casi hasta el final. Hay que meteerlo en el listado de no usarlo en el batch. Escrito en el excel. Ademas estaba mal hecho y lo he tenidoq ue copiar de un backup que tenia yo.
#     - S004_DAY1_MINITWICE_3011_EVENT.csv 234: tiene menos que el DAY2, por el fallo que comentamos.
# - S024_DAY1_MINITWICE_5351: hace fallar al script. Esta es la ultima linea:
#     - S024_MINITWICE_DAY1_5351	889	Quit		9047811	30244` Lo he eliminado y ya esta. Ver abajo. Si anad jiter tb funciona, pero seria falso.
#
# **TRIAD**
# - S022_MINIONCE_5276 > log the Triad Run_I >> Error de Direct X cuando quedaba 1 minuto
# - S005_MINIONCE_2138-vOT_Triads_fMRI_II.log falla, tb tiene un Quit al final: Copiamos las filas del jitter al final y funciona: ademas al final hay que anadir los numeros bien para que la resta de: 10200
# - S007_MINIONCE_4838-vOT_Triads_fMRI_I.log falla por el QUIT en el RUNII tb, por dentro este sujeto es el S011_DAY1_MINITWICE_4838: Copiamos las filas del jitter al final y funciona. Resta ha de dar: 10200
# - Es dificil ajustar: al final voy a elimnar las lineas que pongan quit y ya esta, de hecho, se deberian eleiminar esos volumenes pq no han estado en rest, se les quito la pantalla se supone...
# - El codigo de abajo nos marca correctamente los archivos que no cuadran:
#
# ```
# TRIAD != 283
# S005_MINIONCE_2138_Triad.csv 279: falta solo jitter final
# S007_MINIONCE_4838_Triad.csv 280: falta solo jitter final
# S022_MINIONCE_5276_Triad.csv 268
# S026_DAY2_MINITWICE_2966_Triad.csv 279: falta solo jitter final
# ```
#
#
#
# Ya estan todos los csv creados y revisados para los archivos que tenia hasta el momento. Ahora usando el script de arriba y desde el command line hay que ir convirtiendo los nuevos que lleguen.
# %%
# run first without the if then with the if to see just the ones that are not the required length
# We have to fix the subjects that appear here, it all shuold have the same length for a given test.
# Mark here the ones that we know that have a different length and that we don't care.
print 'BLOCK, print if != 167'
os.chdir(join(fMRI_SPMdir, 'temp_log_files', 'block' ))
A = sorted(glob.glob('*.csv'))
for a in A:
    b = open(a, 'r').readlines()
    if len(b) != 167:
        print a, len(b)
print '\n'

print 'EVENT == variable length, check DAY1 == DAY2'
print 'Print if less than 280 or more than 300'
os.chdir(join(fMRI_SPMdir, 'temp_log_files', 'event' ))
A = sorted(glob.glob('*.csv'))
for a in A:
    b = open(a, 'r').readlines()
    if len(b) <= 280 or len(b) >= 300:
        print a, len(b)
print '\n'

print 'TRIAD, print if != 283'
os.chdir(join(fMRI_SPMdir, 'temp_log_files', 'triad' ))
A = sorted(glob.glob('*.csv'))
for a in A:
    b = open(a, 'r').readlines()
    if len(b) != 283:
        print a, len(b)
print '\n'
print 'END OF CELL'
# %% markdown
# OJO: al mover cada archivo donde le toque habrá que scriptar que cambie el nombre del sujeto por el que le hemos asignado en la bbdd (kepa me ha dicho que lo que aparece dentro del csv no habrá que cambiar, solo el nombre del archivo).
# Crear script python para crear carpetas dentro de data y mover los archivos que toque.
# Copiarlo en vez de moverlo, ya que no ocupa sitio.
#
# Al moverlo cada uno a su sitio, ademas recordar que hay que crear un vector con 1-s de todos los archivos que ya estan para poder marcarlo en el excel.
# %%
# Aqui no hay que cambiar los nombres porque antes de la conversion ya se cambiaron y luego al hacer la conversion
# no es como cuando lee del dicom que vuelve a leer el nombre de la carpeta que esta mal.
GLMdict = dict(zip(GLMs, ['EVENT', 'Triad', 'Block']))
for GLM in GLMs:
    print GLM
    templogglmdir = join(fMRI_SPMdir, 'temp_log_files', GLM)
    os.chdir(templogglmdir)
    todocsv = sorted(glob.glob('*.csv'))
    for csv in todocsv:
        OsirixID = csv.replace('_'+GLMdict[GLM]+'.csv','')
        print OsirixID,csv,
        sub = SEGU['MINID'][SEGU['OsirixID']==OsirixID].tolist()[0] # Cambiar a MINID == sub
        print sub
        source = join(templogglmdir, csv)
        dest = join(fMRI_SPMdir, GLM, 'data', sub, 'log', sub+'.csv')
        if not os.path.isfile(dest):
            # print 'Moving from: ' + source + '\nMoving to: ' + dest + '\n\n'
            shutil.copy2(source, dest)
            shutil.move(source, join('movidos', csv))
        else:
            print csv, 'not copied, exists already'
print '\nEND OF CELL'



# Group Analysis



# Indivdual analysis
# batchPreproc({'S004', 'S005'}, 't', 'event_acpc_lhFusiform', '/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/event')
# cambiar el experiment file y poner arsvct (ya que el paso t crea los c y los ct)
# batchPreproc({'S003','S004', 'S005'}, 'MC', 'event_acpc_lhFusiform', '/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/event')

# Ahora hay que hacer el Analisis
# Primero editar el model_ file para meter la mascara que queramos usar.
# Luego hacer createVEctors, haerlo directamente en Matlab en vez de en Jupyter con qsubs
# createVectors({'S003','S004', 'S005'}, 'event_acpc_lhFusiform')
# Ahora ya podemos lanzar batchAnalyze (ver abajo para lanzarlo con qsubs)
# batchPreproc({'S003','S004', 'S005'}, 'event_acpc_lhFusiform', 'dec', '/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/event')
# %% markdown
# ### Preproc QA
# Hay que ver que este bien alineado
# Sacar fotos y verlo todo junto
#
#
# Crear CMD.txt para cada sujeto con los siguientes datos y lanzar el qsub que los crea todos (igual no se puede pq necesita la pantalla fisicamente, entonces lanzarlos todos serialmente en myCloud)
# -v fMRI_SPM/event/data/S001/anat/highres.nii
# -f freesurferacpc/S001/surf/lh.pial:edgecolor=red
# -f freesurferacpc/S001/surf/lh.white:edgecolor=blue
# -ras -40 -50 -15
# -viewport x
# -ss fMRI_SPM/event/QA/preproc/S001preproc1.png
# -v fMRI_SPM/event/data/S001/func/run1/final/tcvsraf_0222.img:opacity=0.65
# -ss fMRI_SPM/event/QA/preproc/S001preproc2.png
# -quit
# %%
# DO   = True
# SHOW = True
#
# print 'Grabar imagenes para QA'
# GLMs = ['event', 'triad', 'block']
# print 'Obligamos a tener que elegir un GLM cada vez para controlar mejor el output.'
# #######################
# GLM = 'block'
#
# #######################
#
#
# GLMdir = join(fMRI_SPMdir, GLM)
# basedir = ANALYSISdir
#
# for sub in doSubs:
#     OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
#     print '\n\n', sub, OsirixID
#     print '------------------------------'
#
#     cmdText = []
#     cmdText.append('-v ' + join(GLMdir, 'data', sub, 'anat', 'highres.nii'))
#     cmdText.append('-f ' + join(SUBJECTS_DIRacpc, sub, 'surf', 'lh.pial:edgecolor=red'))
#     cmdText.append('-f ' + join(SUBJECTS_DIRacpc, sub, 'surf', 'lh.white:edgecolor=blue '))
#     cmdText.append('-ras -40 -50 -15 ')
#     cmdText.append('-viewport x  ')
#     cmdText.append('-ss ' + join(GLMdir, 'QA', 'preproc', sub + 'Run1preproc1.png'))
#     cmdText.append('-v ' + join(GLMdir, 'data', sub, 'func', 'run1', 'final', 'tcvsraf_0022.img:opacity=0.65'))
#     cmdText.append('-ss ' + join(GLMdir, 'QA', 'preproc', sub + 'Run1preproc2.png'))
#     cmdText.append('-quit')
#
#     subDir = join(GLMdir, 'data', sub)
#     cmdTextFile = open(join(subDir, 'freeviewCmdRun1.txt'), 'w')
#     for text in cmdText:
#         cmdTextFile.write(text + '\n')
#     cmdTextFile.close()
#
#     cmd = str("freeview -cmd " + join(subDir, 'freeviewCmdRun1.txt'))
#     if SHOW: print cmd+'\n'   # Test it before launching
#     if DO: spcmd = sp.call(cmd, shell=True)

# for sub in doSubs:
#     OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
#     print '\n\n', sub, OsirixID
#     print '------------------------------'
#
#     cmdText = []
#     cmdText.append('-v ' + join(GLMdir, 'data', sub, 'anat', 'highres.nii'))
#     cmdText.append('-f ' + join(SUBJECTS_DIRacpc, sub, 'surf', 'lh.pial:edgecolor=red'))
#     cmdText.append('-f ' + join(SUBJECTS_DIRacpc, sub, 'surf', 'lh.white:edgecolor=blue '))
#     cmdText.append('-ras -40 -50 -15 ')
#     cmdText.append('-viewport x  ')
#     cmdText.append('-ss ' + join(GLMdir, 'QA', 'preproc', sub + 'Run2preproc1.png'))
#     cmdText.append('-v ' + join(GLMdir, 'data', sub, 'func', 'run2', 'final', 'tcvsraf_0022.img:opacity=0.65'))
#     cmdText.append('-ss ' + join(GLMdir, 'QA', 'preproc', sub + 'Run2preproc2.png'))
#     cmdText.append('-quit')
#
#     subDir = join(GLMdir, 'data', sub)
#     cmdTextFile = open(join(subDir, 'freeviewCmdRun2.txt'), 'w')
#     for text in cmdText:
#         cmdTextFile.write(text + '\n')
#     cmdTextFile.close()
#
#     cmd = str("freeview -cmd " + join(subDir, 'freeviewCmdRun2.txt'))
#     if SHOW: print cmd+'\n'   # Test it before launching
#     if DO: spcmd = sp.call(cmd, shell=True)
#
# print 'END OF CELL'

# %% markdown
# Ahora hay que crear los vectores, pero eso lo hacemos a mano para cada uno de ellos ya que no tarda nada en matlab.
#
# block todos: hecho
#
# leer todos los sujetos desde groupscurly.m
#
# we have to launch createVectors(all, 'block')
#
# change the model_block for S091 to
#
# MODEL.nscans                =   [184 283]; % 287-4 y dos Runs (), S091= [184 283]
#
# And then launch manually createVectors('S091', 'block') and the batchAnalyze below
#
# qsub -q all.q $mySH/RunMatlab.sh "matlab -nosplash -nodesktop -nodisplay  -r ""batchAnalyze({'S091'},'block','dec','/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/block');exit"""
#
# block S091: hecho
#
#
#
#
# event todos:
# event S044:
# SCANNER.nscans      = [318]; % 322-4 = 318, 1 Run (), S044 = [257]
# createVectors('S044', 'event')
# qsub -q all.q $mySH/RunMatlab.sh "matlab -nosplash -nodesktop -nodisplay  -r ""batchAnalyze({'S044'},'event','dec','/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/event');exit"""
#
#
#
# triad todos: hecho
# triad S022:
# SCANNER.nscans      = [198 198];  % 202-4 = 198, 2 Runs (), S022 = [183 198]
# createVectors('S022', 'triad')
# qsub -q all.q $mySH/RunMatlab.sh "matlab -nosplash -nodesktop -nodisplay  -r ""batchAnalyze({'S022'},'triad','dec','/bcbl/home/public/Gari/MINI/ANALYSIS/fMRI_SPM/triad');exit"""
# %% markdown
# VAmos a convertir todos los spmT a archivos mgh para poder ser visualzados en freesurfer.
# Primero mirar que haya el num de spmT requerido en cada sujeto, con esto (59 para block, 58 para event):
# ```
# find -maxdepth 1 -type d -readable -exec sh -c 'echo "$1"; find "$1"/results/spmT*.img  | wc -l' sh {} ';'
# ```
# Ademas en el mismo paso convertimos toda esta informacion a fsaverage para poder comparar los datos directamente.
# Para esto usare el comando que usa fs ene l recon all para crear el qcache, esto es, tiene los datos de thickness y los pasa a fsaverage, y luego dentro de fsaverage los smoothea con distintos mm. Yo no lo voy a smoothera por ahora hasta ver como queda.
# ```
# mri_surf2surf --srcsubject S001 --srchemi lh --srcsurfreg sphere.reg --trgsubject fsaverage --trghemi lh --trgsurfreg sphere.reg --tval ./tmp.mris_preproc.29656/S001.1.mgh --sval /bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc/S001/surf/lh.thickness --sfmt curv --noreshape --no-cortex
# ```
#
#
#
#

# Moved it from below
# We need to
# %%
DO   = True
SHOW = True

print 'Extract left and right fusiform, and left and right lateral occipital, editado para sumar inferiortemporal'

# Esto me ha dado muchos problemas, en analyze hace un flip que queda mal, con --out_data_type int estropea los numeros, lo dejo en default
GLMs = ['event', 'triad', 'block']
# If you export your dicoms twice over the same subject, all of the dicoms will be duplicated
for sub in doSubs:
    # Find each subject's folder
    OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
    print '\n', sub, OsirixID
    for glm in GLMs:
        aparcFile     = join(SUBJECTS_DIRacpc, sub, 'mri', 'aparc+aseg.mgz')
        aparc2009File = join(SUBJECTS_DIRacpc, sub, 'mri', 'aparc.a2009s+aseg.mgz')

        # cmd1 = str('mri_extract_label '+ aparcFile + ' 1007 ' + join(fMRI_SPMdir, glm, 'data', sub, 'anat', 'lh.fusiform.nii'))
        # qsubcmd1 = str('qsub -v DIR=' + fMRI_SPMdir + ' $mySH/RunPython.sh "' + cmd1 + '"')
        # if SHOW: print qsubcmd1+'\n'   # Test it before launching
        # if DO: spqsubcmd1 = sp.call(qsubcmd1, shell=True)

        # cmd2 = str('mri_extract_label '+ aparcFile + ' 2007 ' + join(fMRI_SPMdir, glm, 'data', sub, 'anat', 'rh.fusiform.nii'))
        # qsubcmd2 = str('qsub -v DIR=' + fMRI_SPMdir + ' $mySH/RunPython.sh "' + cmd2 + '"')
        # if SHOW: print qsubcmd2+'\n'   # Test it before launching
        # if DO: spqsubcmd2 = sp.call(qsubcmd2, shell=True)

        #cmd3 = str('mri_extract_label '+ aparcFile + ' 1008 1031 1029 ' +
        #           join(fMRI_SPMdir, glm,'data', sub, 'anat', 'lh.PPC2.nii'))
        #qsubcmd3 = str('qsub -v DIR=' + SUBJECTS_DIRacpc + ' $mySH/RunPython.sh "' + cmd3 + '"')
        #if SHOW: print qsubcmd3+'\n'   # Test it before launching
        #if DO: spqsubcmd3 = sp.call(qsubcmd3, shell=True)


        #cmd4 = str('mri_extract_label '+ aparc2009File + ' 11159 ' +
        #            join(fMRI_SPMdir, glm, 'data', sub, 'anat', 'lh.ips.nii'))
        #qsubcmd4 = str('qsub -v DIR=' + SUBJECTS_DIRacpc + ' $mySH/RunPython.sh "' + cmd4 + '"')
        #if SHOW: print qsubcmd4+'\n'   # Test it before launching
        #if DO: spqsubcmd4 = sp.call(qsubcmd4, shell=True)


        uno = join(fMRI_SPMdir, glm, 'data', sub, 'anat', 'lh.ips.nii')
        dos = join(fMRI_SPMdir, glm, 'data', sub, 'anat', 'lh.PPC2.nii')
        salida = join(fMRI_SPMdir, glm, 'data', sub, 'anat', 'lh.ipsPPC.nii')
        cmd4 = str('mrcalc ' + uno + ' ' + dos + ' -add - | mrcalc -force - 0 -gt ' + salida)
        qsubcmd4 = str('qsub -v DIR=' + SUBJECTS_DIRacpc + ' $mySH/RunPython.sh "' + cmd4 + '"')
        if SHOW: print qsubcmd4+'\n'   # Test it before launching
        if DO: spqsubcmd4 = sp.call(cmd4, shell=True)
        # if DO: spqsubcmd4 = sp.call(qsubcmd4, shell=True)
print 'END OF CELL'














#
# La idea es luego crear anotaciones y leer informacion de estos archivos tb:
#
# - tamao cluster
# - peak voxel
# - para la annot hay que binarizar y meter en un LUT, segun esto:

# %%
DO   = True
SHOW = True

print 'LANZAR batchAnalyze'
GLMs = ['event', 'triad', 'block']
print 'Obligamos a tener que elegir un GLM cada vez para controlar mejor el output.'
#######################
GLM   = 'block'
modelFilename = 'block_acpc_lhPPC'
steps = 'dec'
#######################
basedir = join(fMRI_SPMdir, GLM)
print doSubs

for sub in doSubs:
    OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
    print '\n\n', sub, OsirixID
    print '------------------------------'
    fcmd = str("batchAnalyze({'" + sub + "'},'" + modelFilename + "','" + steps + "','" +  basedir + "')")
    cmd = str('matlab -nosplash -nodesktop -nodisplay  -r ""' + fcmd + ';exit""')
    qsubcmd = str('qsub -q all.q $mySH/RunMatlab.sh "' + cmd + '"')
    if SHOW: print qsubcmd+'\n'   # Test it before launching
    if not DO: print 'REMEMBER TO CREATE VECTORS'
    if DO: spqsubcmd = sp.call(qsubcmd, shell=True)
print 'END OF CELL'
# %%
# Translate numbers to names (dict copied from $c/MINI/Contrasts.py)
block = dict()
event = dict()
triad = dict()
Translate = dict()
Translate[GLMs[0]] = event
Translate[GLMs[1]] = triad
Translate[GLMs[2]] = block



block['0001'] = 'PSvsNull'
block['0002'] = 'CSvsNull'
block['0003'] = 'PFvsNull'
block['0004'] = 'SDvsNull'
block['0005'] = 'FCvsNull'
block['0006'] = 'PWvsNull'
block['0007'] = 'WLvsNull'
block['0008'] = 'WHvsNull'
block['0009'] = 'CBvsNull'
block['0010'] = 'FFvsNull'
block['0011'] = 'MCTaskvsNull'
block['0012'] = 'endrestvsNull'
block['0013'] = 'ALLvsNull'
block['0014'] = 'AllNoFacesvsNull'
block['0015'] = 'AllFacesvsNull'
block['0016'] = 'RWvsNull'
block['0017'] = 'RWvsPW'
block['0018'] = 'RWvsCS'
block['0019'] = 'RWvsFF'
block['0020'] = 'RWvsSD'
block['0021'] = 'RWvsPS'
block['0022'] = 'RWvsCB'
block['0023'] = 'WHvsPW'
block['0024'] = 'WHvsCS'
block['0025'] = 'WHvsFF'
block['0026'] = 'WHvsSD'
block['0027'] = 'WHvsPS'
block['0028'] = 'WHvsCB'
block['0029'] = 'WLvsPW'
block['0030'] = 'WLvsCS'
block['0031'] = 'WLvsFF'
block['0032'] = 'WLvsSD'
block['0033'] = 'WLvsPS'
block['0034'] = 'WLvsCB'
block['0035'] = 'PWvsCS'
block['0036'] = 'PWvsFF'
block['0037'] = 'PWvsSD'
block['0038'] = 'PWvsPS'
block['0039'] = 'PWvsCB'
block['0040'] = 'CSvsFF'
block['0041'] = 'CSvsSS'
block['0042'] = 'CSvsPS'
block['0043'] = 'CSvsCB'
block['0044'] = 'FFvsSD'
block['0045'] = 'FFvsPS'
block['0046'] = 'FFvsCB'
block['0047'] = 'SDvsPS'
block['0048'] = 'SDvsCB'
block['0049'] = 'PSvsCB'
block['0050'] = 'FacesvsRW'
block['0051'] = 'FacesvsWH'
block['0052'] = 'FacesvsWL'
block['0053'] = 'FacesvsPW'
block['0054'] = 'FacesvsCS'
block['0055'] = 'FacesvsFF'
block['0056'] = 'FacesvsSD'
block['0057'] = 'FacesvsPS'
block['0058'] = 'FacesvsPSFaces'
block['0059'] = 'FacesvsCB'
block['0060'] = 'PWvsRW' # Creado a mano multiplicando por -1 la inversa
block['0061'] = 'RWvsFaces' # Creado a mano multiplicando por -1 la inversa

event['0001'] = 'PSvsNull'
event['0002'] = 'CSvsNull'
event['0003'] = 'PFvsNull'
event['0004'] = 'SDvsNull'
event['0005'] = 'FCvsNull'
event['0006'] = 'PWvsNull'
event['0007'] = 'WLvsNull'
event['0008'] = 'WHvsNull'
event['0009'] = 'CBvsNull'
event['0010'] = 'FFvsNull'
event['0011'] = 'MCTaskvsNull'
event['0012'] = 'ALLvsNull'
event['0013'] = 'AllNoFacevsNull'
event['0014'] = 'AllFacesvsNull'
event['0015'] = 'RWvsNull'
event['0016'] = 'RWvsPW'
event['0017'] = 'RWvsCS'
event['0018'] = 'RWvsFF'
event['0019'] = 'RWvsSD'
event['0020'] = 'RWvsPS'
event['0021'] = 'RWvsCB'
event['0022'] = 'WHvsPW'
event['0023'] = 'WHvsCS'
event['0024'] = 'WHvsFF'
event['0025'] = 'WHvsSD'
event['0026'] = 'WHvsPS'
event['0027'] = 'WHvsCB'
event['0028'] = 'WLvsPW'
event['0029'] = 'WLvsCS'
event['0030'] = 'WLvsFF'
event['0031'] = 'WLvsSD'
event['0032'] = 'WLvsPS'
event['0033'] = 'WLvsCB'
event['0034'] = 'PWvsCS'
event['0035'] = 'PWvsFF'
event['0036'] = 'PWvsSD'
event['0037'] = 'PWvsPS'
event['0038'] = 'PWvsCB'
event['0039'] = 'CSvsFF'
event['0040'] = 'CSvsSD'
event['0041'] = 'CSvsPS'
event['0042'] = 'CSvsCB'
event['0043'] = 'FFvsSD'
event['0044'] = 'FFvsPS'
event['0045'] = 'FFvsCB'
event['0046'] = 'SDvsPS'
event['0047'] = 'SDvsCB'
event['0048'] = 'PSvsCB'
event['0049'] = 'FacesvsRW'
event['0050'] = 'FacesvsWH'
event['0051'] = 'FacesvsWL'
event['0052'] = 'FacesvsPW'
event['0053'] = 'FacesvsCS'
event['0054'] = 'FacesvsFF'
event['0055'] = 'FacesvsSD'
event['0056'] = 'FacesvsPS'
event['0057'] = 'FacesvsPSFaces'
event['0058'] = 'FacesvsCB'
event['0059'] = 'PWvsRW' # Creado a mano multiplicando por -1 la inversa
event['0060'] = 'RWvsFaces' # Creado a mano multiplicando por -1 la inversa


triad['0001'] = 'PWvsNull'
triad['0002'] = 'RWvsNull'
triad['0003'] = 'FFvsNull'
triad['0004'] = 'endrestvsNull'
triad['0005'] = 'ALLvsNull'
triad['0006'] = 'RWPWvsNull'
triad['0007'] = 'AllNonRWvsNull'
triad['0008'] = 'RWvsPW'
triad['0009'] = 'RWvsFF'
triad['0010'] = 'RWvsNonWords'
triad['0011'] = 'PWvsFF'
triad['0012'] = 'SES1_RWvsNonWords'
triad['0013'] = 'SES2_RWvsNonWords'
triad['0014'] = 'ORTH_TO_SES1_RWvsNonWords'
triad['0015'] = 'ORTH_TO_SES2_RWvsNonWords'
triad['0016'] = 'SES1_AllNonRWvsNull'
triad['0017'] = 'SES2_AllNonRWvsNull'
triad['0018'] = 'ORTH_TO_SES1_AllNonRWvsNull'
triad['0019'] = 'ORTH_TO_SES2_AllNonRWvsNull'
# %%
DO   = True
SHOW = True

# Podria crear la inversa del mgh, pero entonces no tendria su equivalente en spmT y en el futuro no
# sabemos para que lo vamos a necesitar. Ahora mismo tal como tengo el codigo no estoy usando el spmT.
print 'Crear contrastes inversos  multiplicando spmT de interes por -1'
GLMs = ['block', 'event']
forBlock = dict()
forEvent = dict()
createContrast = dict()
createContrast['block'] = forBlock
createContrast['event'] = forEvent

# forGLM[newNumber] = ('newContrast', 'oldContrastNumberToInvert')
forBlock['0060']    = ('PWvsRW'   , '0017')
forBlock['0061']    = ('RWvsFaces', '0050')
forEvent['0059']    = ('PWvsRW'   , '0016')
forEvent['0060']    = ('RWvsFaces', '0049')

for GLM in GLMs:
    modelFilename = GLM + '_acpc_lhITfusLatOcc'
    for sub in doSubs:
        OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
        print '\n\n', sub, OsirixID
        print '------------------------------'
        subDir = join(fMRI_SPMdir, GLM, 'analysis_' + modelFilename, 'SUBJECTS', sub)
        os.chdir(join(subDir, 'results'))
        for conNum in createContrast[GLM]:
            cmd1  = str('mrcalc ' + join(subDir, 'results', 'spmT_' + createContrast[GLM][conNum][1] + '.img ') +
                        '-1 -multiply ' +
                        join(subDir, 'results', 'spmT_' + conNum + '.img ')
                        )
            # I am going to to the next to steps in once, this step should dissapear and the new files will be integrated
            cmd2 = str('mri_vol2surf ' +
                       '--srcsubject ' + sub + ' ' +
                       '--projfrac 0.5 ' +
                       '--interp trilinear ' +
                       '--hemi lh ' +
                       '--regheader ' + sub + ' ' +
                       '--mov ' + join(subDir, 'results', 'spmT_' + conNum + '.img') + ' ' +
                       '--o ' + join(subDir, 'results', createContrast[GLM][conNum][0] +'.mgh ')
                       )
            cmd3 = str('mri_surf2surf ' +
                       '--srcsubject ' + sub + ' ' +
                       '--srchemi lh ' +
                       '--srcsurfreg sphere.reg ' +
                       '--sval ' + join(subDir, 'results', createContrast[GLM][conNum][0] +'.mgh') + ' ' +
                       '--trgsubject fsaverage ' +
                       '--trghemi lh ' +
                       '--trgsurfreg sphere.reg ' +
                       '--tval ' + join(subDir, 'results', createContrast[GLM][conNum][0] + '305.mgh') + ' ' +
                       '--sfmt ' +
                       '--curv ' +
                       '--noreshape ' +
                       '--no-cortex '
                       )
            # cmdSuma =  str(cmd1 +  ' && ' +  cmd2 +  ' && '  + cmd3) al qsub no le ha gustado
            qsubcmd = str('qsub -v DIR=' + SUBJECTS_DIRacpc + ' $mySH/RunPython.sh "' + cmd3 + '"')
            if SHOW: print qsubcmd+'\n'   # Test it before launching
            if DO: spqsubcmd = sp.call(qsubcmd, shell=True)

print 'END OF CELL'
# %%
DO   = True
SHOW = True

print 'Crear los MGH desde los spmT'
GLMs = ['event', 'triad', 'block']
print 'Obligamos a tener que elegir un GLM cada vez para controlar mejor el output.'
#######################
GLM   = 'event'
modelFilename = 'event_acpc_lhPPC'
#######################
basedir = join(fMRI_SPMdir, GLM)
print doSubs

for sub in doSubs:
    OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
    print '\n\n', sub, OsirixID
    print '------------------------------'
    subDir = join(fMRI_SPMdir, GLM, 'analysis_' + modelFilename, 'SUBJECTS', sub)
    os.chdir(join(subDir, 'results'))
    spmTresults = glob.glob('spmT*.img')

    for spmT in sorted(spmTresults):
        num = spmT[5:9]
        cmd = str('mri_vol2surf ' +
                  '--srcsubject ' + sub + ' ' +
                  '--projfrac 0.5 ' +
                  '--interp trilinear ' +
                  '--hemi lh ' +
                  '--regheader ' + sub + ' ' +
                  '--mov ' + join(subDir, 'results', spmT) + ' ' +
                  '--o ' + join(subDir, 'results', Translate[GLM][num] +'.mgh')
                 )

        qsubcmd = str('qsub -v DIR=' + SUBJECTS_DIRacpc + ' $mySH/RunPython.sh "' + cmd + '"')
        if SHOW: print qsubcmd+'\n'   # Test it before launching
        if DO: spqsubcmd = sp.call(qsubcmd, shell=True)

print 'END OF CELL'
# %%
DO   = True
SHOW = True

print 'Hacer surf2surf para pasar valores a fsaverage 305, hay que hacerlo en otro paso cuando haya terminado el primero'
GLMs = ['event', 'triad', 'block']
print 'Obligamos a tener que elegir un GLM cada vez para controlar mejor el output.'
#######################
GLM   = 'event'
modelFilename = 'event_acpc_lhPPC'
#######################
basedir = join(fMRI_SPMdir, GLM)
print doSubs

for sub in doSubs:
    OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
    print '\n\n', sub, OsirixID
    print '------------------------------'
    subDir = join(fMRI_SPMdir, GLM, 'analysis_' + modelFilename, 'SUBJECTS', sub)
    os.chdir(join(subDir, 'results'))
    spmTresults = glob.glob('spmT*.img')


    for num in sorted(Translate[GLM]):
        cmd = str('mri_surf2surf ' +
                  '--srcsubject ' + sub + ' ' +
                  '--srchemi lh ' +
                  '--srcsurfreg sphere.reg ' +
                  '--sval ' + join(subDir, 'results', Translate[GLM][num] +'.mgh') + ' ' +
                  '--trgsubject fsaverage ' +
                  '--trghemi lh ' +
                  '--trgsurfreg sphere.reg ' +
                  '--tval ' + join(subDir, 'results', Translate[GLM][num] + '305.mgh') + ' ' +
                  '--sfmt ' +
                  '--curv ' +
                  '--noreshape ' +
                  '--no-cortex'
                 )

        qsubcmd = str('qsub -v DIR=' + SUBJECTS_DIRacpc + ' $mySH/RunPython.sh "' + cmd + '"')
        if SHOW: print qsubcmd+'\n'   # Test it before launching
        if DO: spqsubcmd = sp.call(qsubcmd, shell=True)

print 'END OF CELL'
# %% markdown
# Ahora en matlab obtener las coordenadas de cada contraste.
# Voy a usar el de Jobard para buscar el local.
# Se asume que Jobard es 305, o sea que necesito coordenada local en espacio individuaal:
#
# *****
# I have an RAS from a voxel in MNI305 space (fsaverage space) which I want to convert to MNI152 space. Create a vector of the MNI305 space point like v = [R A S 1]'. Multiply this vector by the matrix below (ie, M*v)
#
#     0.9975   -0.0073    0.0176   -0.0429
#     0.0146    1.0009   -0.0024    1.5496
#    -0.0130   -0.0093    0.9971    1.1840
# Eg, if the RAS point is (10 -20 35), then v = [10 -20 35 1]', and M*v = [10.695 -18.409 36.137 1], so the RAS in MNI152 space would be 10.695 -18.409 36.137.
#
# The above matrix is V152*inv(T152)*R*T305*inv(V305), where V152 and V305 are the vox2ras matrices from the 152 and 305 spaces, T152 and T305 are the tkregister-vox2ras matrices from the 152 and 305 spaces, and R is from $FREESURFER_HOME/average/mni152.register.dat
# *****
#
# O sea que tendre que hacer, para cada sujeto:
# vJobard = [-44 , -58 ,  -15 , 1];
# M = VSubject * inv(TSubject) * R * T305 * inv(V305)
# VSubject =
# *****
# OJO, algo aprecido a esto ya habia hecho en hippovol, pero al reves.
# Le doy la vuelta y ya esta, ver codigo en myGetLocalGlobalMaxima.m
#
# TalXFM = xfm_read([sp filesep 'transforms' filesep 'talairach.xfm']);
# MNI305 = TalXFM * M.vox2ras1 * [i;j;k;1];
#
#
#
# Ahora vamos a crear anotaciones.
# Ya he visto que en freeview puedo ver las superficies, pero no todas juntas:
# ```
# freeview -viewport 3d -f /bcbl/home/public/Gari/MINI/ANALYSIS/freesurferacpc/S006/surf/lh.inflated:annot=aparc.annot:overlay=RWvsCB.mgh:overlay=RWcsNull.mgh:overlay=RWvsPW.mgh:overlay=RWvsFF.mgh:overlay=RWvsCS.mgh:overlay=RWvsSD.mgh:overlay=RWvsPS.mgh &
# ```
# Ademas en la lista Doug explica como crear las anotaciones:
# 1. mri_binarize -- threshold your two raw maps, give one a "--binval 1" the other "--binval 2"
# 2. Add them together: fscalc binmap1.nii add binmap2.nii -o binmap12.nii. This will cause voxels only in binmap1 to be one, only in binmap2 to be 2, and voxels in both to be 3
# 3. Create a lookup table like $FREESURFER_HOME/FreeSurferColorLUT.txt with 4 entries (0, 1, 2, 3), and give them the color you want.
# 4. mris_seg2annot --seg binmap12.nii --ctab yourcolortable ...
#
#  Ver abajo, en vez de cerar una anotiacion que no me sirve tanto, he creado overlays.
#
#
#
#
# ### Analisis de datos
# Ya he conseguido tener toda la informacion en Matlab obtenida en myGetLocalGlobalMaxima y lo tengo en excel. Revisar ese codigo.
# Hay mucha variabilidad. Ahora hay que pensar los siguientes pasos para limpiar estos datos.
# Además, para ir caracterizando mejor los datos, hay que obtener la otra información que aporta fs: girificación, CT, volumen, area, etc. para ver si tenemos algo ahí.
# Para esto tendré que coger el vertex y coger xx vertex contiguos en circulo?
# Creo que lo que haré será lo siguiente:
# 1.- Crear label en el vertex maximo
# 2.- Hacer dilate_label, tener un par de tamaños, ir agrandando, hasta tener un tamaño std reportado en la literatura.
# 3.- Obtener todas las estadisticas en cada uno de los puntos.
# 4.- Además, crear una version de estas mascaras metidas para dentro con el objetivo de ser usadas como ROI en tractografia
# 5.- Tendría sentido obtener la informacion qMRI en: GM, WM pegado a GM, y obviamente en los tractos
#
#
# Ahora tengo que hacer analisis de outliers, etc. Ver si podemos hacer que el GM nunca se vaya más anterior que cierta coordenada...
# PROBLEMAS:
# - casi no hay diferencias en
# - hay una variabilidad brutal.
#
# Voy a crear el label freesurfer
#
# ## UPDATE 25 Oct 2016
# 1. Para reducir la variabilidad habia que hacer algo, ya que tenemos T-s significativas a lo largo de todo el fusiforme, asi qeu mirando Vogel 2012, cogimos el a-c-pVWFA-s (anterior, classical, posterior), y hacer analisis. Me di cuenta de que queaban casi siempre para todos los sujetos fuera del fusiforme, asi que tuve que rehacer los analisis en inferiotemporal, midoccipital y fusiform de lh.aparc.annot de Freesurfer.
# 2. Volvi a hacer analisis R: quitar sujetos con movimiento, quitar outliers, e hice histogramas, scatterplots, medias y SD-s (V06).
# 3. Para poder visualizar los datos cree overlays probabilisticos en fsaverage para cada constraste. Por un lado solo con las maximas y por otro lado creando un ROI con dilate 4 al rededor de la maxima. Los valores de cada overlay son divididos por 97 (num de sujtos) con el objetivo de tener un % en la visualizacion.
# 4. Ahora voy a crear un cmd.txt para sacar una foto a cada uno de estos contrastes. Con esto y la literatura elegiremos los contrastes a usar en cada uno de los ROIs, que se usaran para la busqueda de tractos.

# %%
DO   = True
SHOW = True

print 'Grabar imagenes de los resultados de la max T vertex por cada glm, contraste, ROI'
GLMs = ['event', 'triad', 'block']
print 'Obligamos a tener que elegir un GLM cada vez para controlar mejor el output.'
#######################
GLM = 'block'

#######################


GLMdir = join(fMRI_SPMdir, GLM)
basedir = ANALYSISdir
fsavgDir = join(SUBJECTS_DIRacpc, 'fsaverage')
OverlayDir = join(fsavgDir, 'label', 'MINIprobOverlay')

# Calculate ROI names
dilateLabelBy = ['4', '8', '16', '32']  # Ademas existe el 1 con un vertex
VWFAletter = ['a', 'c', 'p']
autoprefijo = ['GM', 'GMyMin']
for jj in VWFAletter:
    for kk in dilateLabelBy:
        autoprefijo.append(jj + kk)

versionNum = 'V09'
sub = 'fsaverage'
for conName in Translate[GLM].values():
    for prefijo in autoprefijo:
        print '\n\n', GLM, conName, prefijo, versionNum
        print '------------------------------'
        os.chdir(OverlayDir)
        #overlayNameList = glob.glob(GLM+'_'+conName+'_dilBy4_*'+prefijo+versionNum+'.mgh')
        # overlayName = overlayNameList[0]
        overlayName = GLM+'_'+conName+'_'+prefijo+versionNum+'.mgh'

        cmdText = []
        if prefijo[0] == 'G':
            # Iba a poner el mismo thres para todos pero entocnes no se veiran las diferencias
            # Mirare el max across todos y pondre el mismo para todos. Deberia ser 1 pero no se veria nada.
            # El min sera 0.01, que es practicamente 1/100.
            cmdText.append('-f ' + join(SUBJECTS_DIRacpc, sub, 'surf', 'lh.inflated:' +
                                    ':edgethickness=0' +
                                    ':overlay=' + join(OverlayDir, overlayName) +
                                    ':overlay_method=linearopaque'  +
                                    ':overlay_threshold=0.01,0.8 '
                                    ))
        else:
            labelName = prefijo[0] + 'VWFA' + prefijo[1:len(prefijo)] + '.label'
            cmdText.append('-f ' + join(SUBJECTS_DIRacpc, sub, 'surf', 'lh.inflated:' +
                                    ':edgethickness=0' +
                                    ':label=' + join(fsavgDir, 'label', labelName) +
                                    ':label_outline=true'
                                    ':overlay=' + join(OverlayDir, overlayName) +
                                    ':overlay_method=linearopaque' +
                                    ':overlay_threshold=0.01,0.8 '
                                    ))

        cmdText.append('-viewport 3D  ')
        cmdText.append('-cam dolly 1.5 elevation -60 azimuth 35 ')
        cmdText.append('-ss ' + join(GLMdir, 'QA', 'MaxVertexProbMap',
                       'ProbMap_'+GLM+'_'+conName+'_'+prefijo+versionNum+'.png'))
        cmdText.append('-quit')

        subDir = join(GLMdir, 'data', sub)
        cmdTextFile = open(join(OverlayDir, 'cmd_'+GLM+'_'+conName+'_'+prefijo+versionNum+'.txt'), 'w')
        for text in cmdText:
            cmdTextFile.write(text + '\n')
        cmdTextFile.close()

        cmd = str("freeview -cmd " + join(OverlayDir, 'cmd_'+GLM+'_'+conName+'_'+prefijo+versionNum+'.txt'))
        if SHOW: print cmd+'\n'   # Test it before launching
        if DO: spcmd = sp.call(cmd, shell=True)

print 'END OF CELL'
# %% markdown
# # TRIADS
# Check in the Notes app the fotographs of the Seghier paper analysis.
# The summary is that we want the SM network activated (WORDS > NonWords) in the DefaultNetwork (NULL > NonWords) in the Angular area.
# In the paper he defines a little bit differently the areas, but the most ventral zone in which we are interested has both SM > Fixation and Fixation > PM, being PM the perceptual network. In our case the PM is NonWord > NULL, and precisely in our case this is the Default Network.
#
#
# In order to do the single subject analysis we need to specify a previous fROI and if we want a mask as well.
# The idea is going to be:
# 1. Create a square ROI mask around the local maximas. The important thing here is that as we can see below, the overlap center reported by Seghier is the same of our local/global maxima (depending on how we did it, but it is basically the same). So, we are going to create a square around this local maxima with the idea of limiting the analyses to that square.
#
# SQUARE MASK FOR TRIADS:
# center for Seghier: -48        -68                28
# Center for us: -45.5/-48      -69.5/-67        25/27.5
# Data for the square: x = ; y = ; z =
#
# 2. Run again the spm_ss toolbox to check
#     1. Effect(s) of interest: Word-NonWord
#     2. Localizer: NonWord-Null
#     3. Con esto lo que pasa es que no son orthogonales, pq NonWords estan en todas, por lo que spm_ss pega un aviso y luego
#
# 3.- SPM_SS
#     1. Buscar cual es el local maxima para cada sujeto.
#     2. El xyz = [-48, -67,  27.5] es el del general.
#     3. Buscar funcion en xjview.m o SPM que nos haga esto.
#     4. Asi ponemos lo que nos interesa: xyz=spm_mip_ui('SetCoords',xyz,hMIPax)
#     5.
#
#
#
#
#     DefNet: pValue = 0.001, parte z min = 17, x min = -38
#     Crear caja para enmascarar el Def Network: x: -25 -90, y: -50 -120, z: 10 100
#
#
# Hay que renombrar analysis_triad a analysis_triad_OLD, y luego volver a crear todos los analisis con el mask.
# 1.- Hacer de nuevo el createVectors (ojo, cada modelo tiene un sujeto que hay que hacer aparte editando experiment.m).
# 2.- Lanzar todos los qsubs con el script de arriba para lanzar el 1st level analisis.

# %% markdown
# ## OTROS
# En el proceso de analysis hay que ir haciendo tareas para todos los sujetos, ir poniendo aqui para ir aprovechando los codigos.
# %%
# Remove duplicated subjects dicomsfor sub in doSubs:
    # Find each subject's folder
    OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
    print '\n', sub, OsirixID
# If you export your dicoms twice over the same subject, all of the dicoms will be duplicated
for sub in doSubs:
    # Find each subject's folder
    OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
    print '\n', sub, OsirixID
    i = 0
    os.chdir(join(dicomdir, OsirixID))
    folders = glob.glob('*')
    for a in folders:
        os.chdir(join(dicomdir, OsirixID, a))
        tipos = glob.glob('*001.dcm')
        cuantos = len(tipos)
        # if cuantos >= 2:
        if cuantos == 2:
            #for tipo  in tipos:
                #longitud = len(tipo.split('-'))
                #if longitud == 3:
            guardarGrep = 'IM-'+tipos[0].split('-')[1]
            borrarGrep = 'IM-'+tipos[1].split('-')[1]
            i += 1
            print i, join(dicomdir, OsirixID, a, tipos[0]), guardarGrep
            print join(dicomdir, OsirixID, a, tipos[1]), borrarGrep
        else:
            multi = glob.glob('*-0001-0001.dcm')
            if len(multi) == 2:
                guardarGrep = 'IM-'+multi[0].split('-')[1]
                borrarGrep = 'IM-'+multi[1].split('-')[1]
                i += 1
                print i, join(dicomdir, OsirixID, a, multi[0]), guardarGrep
                print join(dicomdir, OsirixID, a, multi[1]), borrarGrep
        for f in os.listdir('./'):
            if re.search(borrarGrep, f):
                print 'REMOVE: ', join(dicomdir, OsirixID,a,f)
                #os.remove(join(dicomdir, OsirixID,a,f))
print 'END OF CELL'
# %%
DO = True
SHOW = True

print 'Rename all anat folder names to anat_noAcpc'
GLMs = ['event', 'triad', 'block']
# If you export your dicoms twice over the same subject, all of the dicoms will be duplicated
for sub in doSubs:
    # Find each subject's folder
    OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
    print '\n', sub, OsirixID
    for glm in GLMs:
        srcDir = join(fMRI_SPMdir, glm, 'data', sub, 'anat')
        dstDir = join(fMRI_SPMdir, glm, 'data', sub, 'anat_noAcpc')
        if SHOW: print srcDir, dstDir, '\n'
        if DO: shutil.move(srcDir, dstDir)

print 'END OF CELL'

# %%
DO = True
SHOW = True

print 'Copy rawavg.mgz to nii format to the event/subName/data/anat'

# Esto me ha dado muchos problemas, en analyze hace un flip que queda mal, con --out_data_type int estropea los numeros, lo dejo en default
GLMs = ['event', 'triad', 'block']
# If you export your dicoms twice over the same subject, all of the dicoms will be duplicated
for sub in doSubs:
    # Find each subject's folder
    OsirixID = segu['OsirixID'][segu['MINID']==sub].tolist()[0] # Cambiar a OsirixID
    print '\n', sub, OsirixID
    for glm in GLMs:
        srcFile   = join(SUBJECTS_DIRacpc, sub, 'mri', 'rawavg.mgz')
        dstFile   = join(fMRI_SPMdir, glm, 'data', sub, 'anat', 'highres.nii')

        # if DO: os.mkdir(join(fMRI_SPMdir, glm, 'data', sub, 'anat'))

        # cmd = str('mri_convert -i '+ srcFile + ' -o ' + dstFile + ' --out_data_type int')
        cmd = str('mri_convert -i '+ srcFile + ' -o ' + dstFile)
        qsubcmd = str('qsub -v DIR=' + fMRI_SPMdir + ' $mySH/RunPython.sh "' + cmd + '"')
        if SHOW: print qsubcmd+'\n'   # Test it before launching
        if DO: spqsubcmd = sp.call(qsubcmd, shell=True)

print 'END OF CELL'
# %%
DO   = False
SHOW = True

print 'Grabar imagenes de los resultados de la max T vertex por cada glm, contraste, ROI'
GLMs = ['event', 'triad', 'block']
print 'Obligamos a tener que elegir un GLM cada vez para controlar mejor el output.'
#######################
GLM = 'block'

#######################


GLMdir = join(fMRI_SPMdir, GLM)
basedir = ANALYSISdir
fsavgDir = join(SUBJECTS_DIRacpc, 'fsaverage')
OverlayDir = join(fsavgDir, 'label', 'MINIprobOverlay')

# Calculate ROI names
dilateLabelBy = ['4', '8', '16', '32']  # Ademas existe el 1 con un vertex
VWFAletter = ['a', 'c', 'p']
autoprefijo = ['GM', 'GMyMin']
for jj in VWFAletter:
    for kk in dilateLabelBy:
        autoprefijo.append(jj + kk)

versionNum = 'V09'
sub = 'fsaverage'
for conName in Translate[GLM].values():
    for prefijo in autoprefijo:
        print '\n\n', GLM, conName, prefijo, versionNum
        print '------------------------------'
        os.chdir(OverlayDir)
        #overlayNameList = glob.glob(GLM+'_'+conName+'_dilBy4_*'+prefijo+versionNum+'.mgh')
        # overlayName = overlayNameList[0]
        overlayName = GLM+'_'+conName+'_'+prefijo+versionNum+'.mgh'

        cmdText = []
        if prefijo[0] == 'G':
            # Iba a poner el mismo thres para todos pero entocnes no se veiran las diferencias
            # Mirare el max across todos y pondre el mismo para todos. Deberia ser 1 pero no se veria nada.
            # El min sera 0.01, que es practicamente 1/100.
            cmdText.append('-f ' + join(SUBJECTS_DIRacpc, sub, 'surf', 'lh.inflated:' +
                                    ':edgethickness=0' +
                                    ':overlay=' + join(OverlayDir, overlayName) +
                                    ':overlay_method=linearopaque'  +
                                    ':overlay_threshold=0.01,0.8 '
                                    ))
        else:
            labelName = prefijo[0] + 'VWFA' + prefijo[1:len(prefijo)] + '.label'
            cmdText.append('-f ' + join(SUBJECTS_DIRacpc, sub, 'surf', 'lh.inflated:' +
                                    ':edgethickness=0' +
                                    ':label=' + join(fsavgDir, 'label', labelName) +
                                    ':label_outline=true'
                                    ':overlay=' + join(OverlayDir, overlayName) +
                                    ':overlay_method=linearopaque' +
                                    ':overlay_threshold=0.01,0.8 '
                                    ))

        cmdText.append('-viewport 3D  ')
        cmdText.append('-cam dolly 1.5 elevation -60 azimuth 35 ')
        cmdText.append('-ss ' + join(GLMdir, 'QA', 'MaxVertexProbMap',
                       'ProbMap_'+GLM+'_'+conName+'_'+prefijo+versionNum+'.png'))
        cmdText.append('-quit')

        subDir = join(GLMdir, 'data', sub)
        cmdTextFile = open(join(OverlayDir, 'cmd_'+GLM+'_'+conName+'_'+prefijo+versionNum+'.txt'), 'w')
        for text in cmdText:
            cmdTextFile.write(text + '\n')
        cmdTextFile.close()

        cmd = str("freeview -cmd " + join(OverlayDir, 'cmd_'+GLM+'_'+conName+'_'+prefijo+versionNum+'.txt'))
        if SHOW: print cmd+'\n'   # Test it before launching
        if DO: spcmd = sp.call(cmd, shell=True)

print 'END OF CELL'
