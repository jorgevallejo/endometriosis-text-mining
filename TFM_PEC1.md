# Trabajo de Fin de Máster - PEC1 Plan de Trabajo  \
## Jorge Vallejo Ortega
## 08/03/2021

# Endo-Mining: herramienta web para la búsqueda automatizada de
genes potencialmente relacionados con la endometriosis a través de
minería de textos en PubMed

##1. ​ Contexto y Justificación del Trabajo

###1.1. Descripción general:

Se desarrollará una aplicación web interactiva para descubrir información contenida en los abstracts de las publicaciones científicas en Pubmed mediante el uso de minería de textos. Se aplicará a la búsqueda de genes relacionados con la endometriosis.

###1.2 Justificación del TFG

La endometriosis es una enfermedad del sistema reproductor femenino caracterizada por crecimiento de tejido endometrial ectópico fuera del útero. Afecta a alrededor del 10% de mujeres en edad reproductora (Malvezzi et al. 2020). Sus principales síntomas son infertilidad y dolor (Rolla 2019), y se ha asociado con depresión y fatiga(Chapron et al., 2019); por lo que supone una importante disminución en su calidad de vida del paciente.

No existe ninguna prueba genética que identifique personas con mayor riesgo de desarrollar la enfermedad, ni biomarcadores que permitan un diagnóstico fiable. Sin embargo podemos acceder a una amplia literatura científica que explora la relación entre genes y endometriosis (Rolla 2019) y, aunque se estima una heredabilidad de entre 0.27 a 0.51, no se ha desmostrado relación directa entre este y ninguno de los genes candidatos. Se le supone una etiología multifactorial en la que intervendrían diferentes factores genéticos y ambientales (Sapkota, 2017).

PubMed es un recurso en línea de acceso público y gratuito consistente en una base de datos - en continuo crecimiento - que incluye más de 32 millones de citas y abstracts de literatura biomédica, tanto artículos (MEDLINE y PubMed Central) como libros (Bookshelf)(dos citas web de la National Library of Medicine). La inmensa cantidad de información contenida en la base de datos hace que incluso búsquedas restrictivas recuperen muchas veces cientos (o miles) de citas relevantes, un flujo de información difícilmente digerible por el investigador mediante los métodos tradicionales de lectura y análisis de artículos individuales.

Las técnicas de minería de textos permiten a los investigadores de las áreas biomédicas un acceso efectivo y eficiente al conocimiento enterrado en las ingentes cantidades de literatura publicada, además de suplementar la información extraída mediante minería de datos a partir de otras fuentes de datos masivos como la secuenciación de genomas, datos de expresión génica y de estructuras proteicas. Ambas funciones permiten acelerar la investigación biomédica (Aggarwal, 2012).

En resumen, en este trabajo estamos abordando el estudio de una enfermedad - la endometriosis - que afecta negativamente la calidad de vida de un importante porcentaje de la población mundial. Lo hacemos centrándonos en un aspecto todavía muy desconocido como es su relación con la genética, y usando una herramienta analítica y exploratoria - la minería de textos - que permite extraer información de una enorme fuente de datos de acceso público pero poco estructurada como son los abstracts de los artículos científicos. Desarrollando una aplicación web que permita a cualquier persona reproducir el proceso de forma sencilla y automática, estamos contribuyendo dentro de nuestras posibilidades al crecimiento del conocimiento científico y técnico al exponer de forma pública los resultados de nuestro trabajo y el modo en el que hemos llegado hasta ellos.

##2. Objetivos

###2.1 Objetivo general

* Encontrar genes relacionados con la endometriosis aplicando técnicas de minería de textos.

###2.2 Objetivos específicos

1. Desarrollar un script que permita realizar un procedimiento de minería de textos automáticamente, desde la recopilación de datos en bruto hasta la presentación de resultados.

2. Desarrollar una aplicación web implementando el script de minería de textos resultante del objetivo anterior.

##3. Enfoque y método a seguir

La herramienta se desarrollará utilizando el lenguaje de programación R. En primer lugar, debido a que dispone una amplia variedad de librerías especialmente enfocadas a la minería de textos en general (lsa, tidytext, tm) y al acceso a datos de PubMed y el NCBI (easyPubmed, pubmed.mineR, bibliometrix, rentrez, RISmed). Y, en segundo lugar, debido a que es el lenguaje mejor conocido por el escritor de este trabajo.

El flujo de trabajo se basará en el expuesto en Rani _et al._ (2015) y Liu (2016), y constará de los siguientes pasos:

1. Reducción de la información. En esta fase se descargará un pequeño subgrupo de abstracts de todos los disponibles en la base de datos PubMed. Para seleccionar dichos abstracts de interés se usarán palabras clave apropiadas al tema de este trabajo y se acotará por fechas. Estos abstracts conformarán el corpus primario de documentos sobre los que se realizará la minería de textos.

2. Preprocesado. Durante esta fase procederemos a la itemización de los textos que componen el corpus primario. Ésta consistirá en el desglose de cada abstract en las frases y palabras que lo componen.

3. Reconocimiento y normalización de entidades. Identificaremos aquellos ítems correspondientes a nombres de genes; consolidaremos además sinónimos de genes a un identificador único de cada gen. El resultado de este paso es una lista de genes posiblemente asociados con la endometriosis según la literatura biomédica que compone el corpus primario.

4. Filtrado estadístico. Se realizará un filtrado de la lista mediante un test frente a la distribución hipergeométrica, para distinguir aquellos genes con menor probabilidad de haber sido recuperados por azar.

5. Caracterización de la lista final de genes. La lista de genes candidatos se examinará buscando categorías funcionales enriquecidas según los términos de ontología génica (GO, _gene ontology_ ). Se tendrán en cuenta las categorías de proceso biológico, componente celular y función molecular.

6. Visualización de la información en forma de tablas y gráficas (nubes de palabras con los genes, diagramas de barras de términos GO).

Finalmente, el flujo de trabajo descrito se integrará en el desarrollo de la aplicación web usando el paquete Shiny de R.