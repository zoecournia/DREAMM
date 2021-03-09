<?xml version="1.0" encoding="ISO-8859-1"?>  
<helpset>  
    <title>Descriptor Specifications</title>  
    <maps>  
        <!-- Página por defecto al mostrar la ayuda -->  
        <homeID>gcf</homeID>  
        <!-- Que mapa deseamos -->  
        <mapref location="urlsidx.jhm" />  
    </maps>  
  
    <!-- Las Vistas que deseamos mostrar en la ayuda -->  
    <view>  
        <!-- Deseamos una tabla de contenidos -->  
        <name>Content</name>  
        <!-- El tooltiptext de la vista -->  
        <label>Content Table </label> 				
        <type>javax.help.TOCView</type>  
        <!-- El icono que se muesta -->  
        <image></image>  
        <!-- El fichero que la define -->  
        <data>treeidx.xml</data>  
    </view>  
  
    <view xml:lang="en">  
        <!-- Deseamos que se puedan realizar búsquedas -->  
        <name>Search</name>  
        <!-- El tooltiptext -->  
        <label>Search</label>  
        <!-- El icono que se muesta -->  
        <image>find</image>  
        <type>javax.help.SearchView</type>  
        <data engine="com.sun.java.help.search.DefaultSearchEngine">  
            JavaHelpSearch  
        </data>  
    </view>  
  
    <!-- Definición de la ventana principal de la ayuda-->  
    <presentation default="true" displayviews="true" displayviewimages="true">  
        <name>MainWin</name>  
        <!-- Dimensiones iniciales -->  
        <size width="640" height="480" />   
        <!-- Posición inicial -->  
        <location x="200" y="200" />  
        <!-- Título de la ventana -->  
        <title>Descriptor Theory</title> 
        <!-- icono de la ventana -->  
        <image>window</image>
        <!-- Definimos la barra de herramientas de la ventana -->  
        <toolbar>  
            <!-- Permitimos ir a la página anterior -->  
            <helpaction image="previuos" >  
                javax.help.BackAction  
            </helpaction>  
            <!-- Permitimos ir a la página siguiente -->  
            <helpaction image="next" >  
                javax.help.ForwardAction  
            </helpaction>  
            <!-- Permitimos imprimir el contenido   
            <helpaction image="print" >  
                javax.help.PrintAction  
            </helpaction>  
            Permitimos configurar la impresión   
            <helpaction image="printsetup" >  
                javax.help.PrintSetupAction  
            </helpaction>  -->
        </toolbar>  
    </presentation>  
</helpset>  