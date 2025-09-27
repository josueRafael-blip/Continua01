import streamlit as st
import itertools
import os
import zipfile
import tempfile
import io
import sys

# Configuración para evitar warnings de RDKit
import warnings
warnings.filterwarnings('ignore')

# Manejo de importación de RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    st.error("❌ RDKit no está instalado. Por favor instala RDKit para usar la funcionalidad de conversión a XYZ.")
    st.info("Instala con: pip install rdkit")
    RDKIT_AVAILABLE = False

# --- (MANTENER TODAS TUS FUNCIONES DE QUÍMICA SIN CAMBIOS) ---
def detectar_quiralidad(smiles: str):
    if not RDKIT_AVAILABLE:
        return False, "RDKit no disponible", []
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "SMILES inválido", []
            
        centros = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            
        if len(centros) == 0:
            return False, "Su molécula no es quiral", []
        else:
            return True, f"Su molécula es quiral. Se detectaron {len(centros)} posibles centros", centros
             
    except Exception as e:
        return False, f"Error al analizar la molécula: {str(e)}", []

def analizar_centros_existentes(smiles: str):
    centros_especificados = 0
    posiciones_at = []
    i = 0
    
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                centros_especificados += 1
                posiciones_at.append(i)
                i += 2
            else:
                centros_especificados += 1
                posiciones_at.append(i)
                i += 1
        else:
            i += 1
            
    return centros_especificados, posiciones_at

def generar_estereoisomeros(smiles: str):
    posiciones = []
    i = 0
    while i < len(smiles):
        if smiles[i] == "@":
            if i + 1 < len(smiles) and smiles[i+1] == "@":
                posiciones.append((i, True))
                i += 2
            else:
                posiciones.append((i, False))
                i += 1
        else:
            i += 1
            
    n = len(posiciones)
    
    if n == 0:
        st.warning("⚠️ El SMILES no tiene centros quirales especificados con @ o @@. No se generarán isómeros.")
        return [], n
    elif n > 3:
        st.error("❌ El SMILES tiene más de 3 centros quirales. No se generarán isómeros.")
        return [], n
    
    combinaciones = list(itertools.product(["@", "@@"], repeat=n))
    resultados = []
    
    for comb in combinaciones:
        chars = list(smiles)
        offset = 0
        for (pos, era_doble), val in zip(posiciones, comb):
            real_pos = pos + offset
            if era_doble:
                chars[real_pos:real_pos+2] = list(val)
                offset += len(val) - 2
            else:
                chars[real_pos:real_pos+1] = list(val)
                offset += len(val) - 1
        resultados.append("".join(chars))
            
    return resultados, n


def smiles_to_xyz(smiles, mol_id):
    if not RDKIT_AVAILABLE:
        return None, "❌ RDKit no está disponible"
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, f"❌ Error: SMILES inválido {smiles}"
            
        mol = Chem.AddHs(mol)
            
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  
            
        embed_result = AllChem.EmbedMolecule(mol, params)
        if embed_result != 0:
            params.useRandomCoords = True
            embed_result = AllChem.EmbedMolecule(mol, params)
            if embed_result != 0:
                return None, f"⚠️ No se pudo generar conformación 3D para {smiles}"
            
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except:
            pass
            
        conf = mol.GetConformer()
        xyz_content = f"{mol.GetNumAtoms()}\n{smiles}\n"
            
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            xyz_content += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
            
        return xyz_content, f"✅ Molécula {mol_id} procesada correctamente"
            
    except Exception as e:
        return None, f"❌ Error procesando {smiles}: {str(e)}"

def crear_archivo_zip(archivos_xyz):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for filename, content in archivos_xyz.items():
            zip_file.writestr(filename, content)
    return zip_buffer.getvalue()
# -----------------------------------------------------


def main():
    # --- CONFIGURACIÓN DE PÁGINA ---
    st.set_page_config(
        page_title="Inchiral - Generador de Estereoisómeros",
        page_icon="🧪",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    # --- NUEVO ESTILO: TEMA OSCURO 'LABORATORIO AVANZADO' ---
    PRIMARY_COLOR = "#00bfa5"  # Teal brillante/Aqua
    SECONDARY_COLOR = "#9c27b0" # Violeta
    BACKGROUND_COLOR = "#212121" # Gris oscuro casi negro
    TEXT_COLOR = "#e0e0e0" # Gris claro

    st.markdown(f"""
        <style>
        /* Estilo Principal (Fondo oscuro) */
        .main {{
            background-color: {BACKGROUND_COLOR};
            color: {TEXT_COLOR};
        }}
        /* Sidebar (Fondo un poco más claro) */
        .css-1d3f82m, .st-emotion-cache-1c5y2fy {{ 
            background-color: #333333; 
            color: {TEXT_COLOR};
            border-right: 2px solid {PRIMARY_COLOR};
        }}
        
        /* Títulos */
        h1, h2, h3, h4, .st-emotion-cache-10tr2k1.e1ezn98l1 {{
            color: {PRIMARY_COLOR}; 
            font-family: 'Consolas', 'Courier New', monospace; /* Fuente técnica/científica */
        }}
        
        /* Texto General */
        .stMarkdown, .stText, .stTextInput > label, .stSelectbox > label {{
            color: {TEXT_COLOR};
        }}

        /* Botones (Primario con color vibrante) */
        .stButton>button, .stDownloadButton>button {{
            background-color: {SECONDARY_COLOR}; 
            color: white;
            border-radius: 6px;
            border: 1px solid {SECONDARY_COLOR};
            padding: 10px 20px;
            font-weight: bold;
            transition: background-color 0.3s;
        }}
        .stButton>button:hover, .stDownloadButton>button:hover {{
            background-color: #7b1fa2; /* Tono más oscuro de violeta */
            color: white;
        }}
        
        /* Campos de Entrada de Texto */
        .stTextInput>div>div>input {{
            border-radius: 6px;
            border: 1px solid #666666;
            background-color: #383838; /* Fondo oscuro en input */
            color: {TEXT_COLOR};
            padding: 10px;
        }}
        
        /* Mensajes de Estado (Colores de neón/futuristas) */
        .stAlert {{
            border-radius: 6px;
            padding: 15px;
        }}
        .stAlert.st-info {{
            background-color: #004d40; /* Fondo Teal Oscuro */
            color: {PRIMARY_COLOR};
            border-left: 5px solid {PRIMARY_COLOR};
        }}
        .stAlert.st-success {{
            background-color: #1b5e20; /* Fondo Verde Oscuro */
            color: #69f0ae; /* Verde Neón */
            border-left: 5px solid #69f0ae;
        }}
        .stAlert.st-warning {{
            background-color: #5d4037; /* Fondo Naranja/Marrón Oscuro */
            color: #ffb74d; /* Naranja Brillante */
            border-left: 5px solid #ffb74d;
        }}
        .stAlert.st-error {{
            background-color: #b71c1c; /* Fondo Rojo Oscuro */
            color: #ff8a80; /* Rojo Claro */
            border-left: 5px solid #ff8a80;
        }}
        
        /* Código (SMILES) */
        .stCode {{
            background-color: #383838; /* Gris oscuro en código */
            border: 1px solid #555555;
            border-left: 5px solid {SECONDARY_COLOR};
            color: #ffe082; /* Color amarillo/cálido para el código */
            border-radius: 5px;
        }}
        
        /* Contenedores y Expander */
        .stContainer, .stExpander {{
            border: 1px solid #383838;
            border-radius: 10px;
            padding: 15px;
            margin-bottom: 15px;
            background-color: #2c2c2c; /* Fondo de contenedores */
        }}
        </style>
        """, unsafe_allow_html=True)

    # -----------------------------------------------------

    st.title("🧪 **Inchiral** | Generador de Estereoisómeros")
    st.markdown(f'<p style="font-size: 18px; color: {PRIMARY_COLOR};"><strong>Análisis y manipulación avanzada de quiralidad molecular.</strong></p>', unsafe_allow_html=True)
    st.markdown("---")
    
    # 3. Sidebar (Ajustado al nuevo tema)
    with st.sidebar:
        try:
            # Puedes intentar una imagen monocromática para que encaje mejor con el tema oscuro
            st.image("imagenes1/logo.png", width=200) 
        except:
            st.markdown(f"## 🔬 **Inchiral**", unsafe_allow_html=True)
            
        st.markdown("---")
        st.subheader("🛠️ **Instrucciones**")
        st.markdown(f'<p style="color: {TEXT_COLOR}">Este módulo utiliza RDKit para detectar centros quirales y generar todas las combinaciones estéreo-isoméricas (máx. 3 centros).</p>', unsafe_allow_html=True)
        st.markdown("""
        **Flujo de trabajo Y Ejemplos:**
        1. Ingresa un código **SMILES** (puede incluir o no quiralidad).  
           - Ejemplo sin centros quirales: `C=CC`  
           - Ejemplo con un centro quiral: `CC(O)C`  
        2. El sistema detecta automáticamente los **centros quirales**.  
           - Ejemplo con estereoquímica definida: `C[C@H](N)O`  
        3. Si hay estereoquímica definida con `@` o `@@`, se generan los posibles estereoisómeros.  
           - Ejemplo de aminoácido (alanina): `N[C@H](C)C(=O)O`  

        """)
        st.markdown("---")
        st.markdown('<div style="text-align: center; color: #555555;"><small>Interfaz de Laboratorio v2.0</small></div>', unsafe_allow_html=True)
    
    # 4. Contenido Principal
    st.subheader("📝 Entrada de Datos SMILES")
    
    with st.container():
        smiles_input = st.text_input(
            "👉 Ingresa el código SMILES:",
            placeholder="Ejemplo: C[C@H](O)[C@@H](N)C"
        )
        
    st.markdown("---")

    if smiles_input:
        st.subheader("🔬 Resultados del Análisis Quiral")
        
        es_quiral, mensaje_quiralidad, centros_detectados = detectar_quiralidad(smiles_input)
        centros_especificados, posiciones_at = analizar_centros_existentes(smiles_input)

        col1, col2 = st.columns(2)
        
        with col1:
            with st.container():
                st.markdown(f"**🔎 Análisis de Molécula (RDKit):**")
                if RDKIT_AVAILABLE:
                    if es_quiral:
                        st.success(f"✅ {mensaje_quiralidad}")
                        if centros_detectados:
                            st.markdown("• **Índices de átomos quirales:**")
                            for i, (idx, chirality) in enumerate(centros_detectados):
                                st.markdown(f"&nbsp;&nbsp;&nbsp;&nbsp;• Átomo **{idx}**")
                    else:
                        st.warning(f"⚠️ {mensaje_quiralidad}")
                else:
                    st.error("❌ RDKit no disponible")
        
        with col2:
            with st.container():
                st.markdown(f"**📋 Centros para Generación Estérica:**")
                if centros_especificados > 0:
                    st.info(f"✅ **{centros_especificados}** centros definidos (@/@ @). Generando **{2**centros_especificados}** isómeros.")
                    st.markdown("• **Posiciones de marcadores @:**")
                    st.code(str(posiciones_at))
                else:
                    st.warning("⚠️ No hay centros especificados. No se generarán isómeros.")

        st.markdown("---")

        if RDKIT_AVAILABLE and es_quiral and centros_especificados == 0:
            st.info("""
            **💡 Alerta:** La molécula es quiral, pero **no** tiene la estereoquímica definida en el SMILES.
            Añade **@** o **@ @** a los centros para generar los isómeros.
            """)
            st.markdown("---")

        isomeros, n_centros = [], 0
        if centros_especificados > 0:
            with st.spinner(f"🔄 **Protocolo Activo:** Generando {2**centros_especificados} estereoisómeros..."):
                isomeros, n_centros = generar_estereoisomeros(smiles_input)
        
        tab1, tab2, tab3 = st.tabs(["📋 Lista Completa (SMILES)", "💾 Exportar .smi", "⚛️ Generar Coordenadas XYZ"])
        
        if isomeros:
            # TAB 1: Lista
            with tab1:
                st.markdown(f"### **Total:** {len(isomeros)} Isómeros Estéricos")
                colA, colB = st.columns(2)
                for i, isomero in enumerate(isomeros):
                    col = colA if i % 2 == 0 else colB
                    col.code(f"MOL-{i+1}. {isomero}")
            
            # TAB 2: Descargar SMI
            with tab2:
                smi_content = "\n".join(isomeros)
                st.download_button(
                    label="📥 Descargar Archivo .smi",
                    data=smi_content,
                    file_name="estereoisomeros.smi",
                    mime="text/plain",
                    key="download_smi_button"
                )
                with st.expander("👀 Ver Archivo de Texto"):
                    st.text(smi_content)
            
            # TAB 3: Convertir a XYZ
            with tab3:
                if st.button("🚀 Iniciar Conversión a XYZ (Generación 3D)", type="primary", key="convert_xyz_button"):
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    archivos_xyz = {}
                    mensajes_log = []
                    
                    for i, smiles in enumerate(isomeros):
                        try:
                            progress = (i + 1) / len(isomeros)
                            progress_bar.progress(progress)
                            status_text.text(f"Procesando MOL-{i+1}/{len(isomeros)}: {smiles}")
                            
                            xyz_content, mensaje = smiles_to_xyz(smiles, i+1)
                            mensajes_log.append(mensaje)
                            
                            if xyz_content:
                                archivos_xyz[f"mol_{i+1}.xyz"] = xyz_content
                        except Exception as e:
                            mensajes_log.append(f"❌ Error procesando MOL-{i+1}: {str(e)}")
                            
                    progress_bar.progress(1.0)
                    status_text.success("✅ **Protocolo Finalizado:** Archivos 3D generados.")
                    
                    with st.expander("📋 Log Detallado de Conversión"):
                        for mensaje in mensajes_log:
                            if "❌" in mensaje or "⚠️" in mensaje:
                                st.error(mensaje)
                            else:
                                st.success(mensaje)
                    
                    if archivos_xyz:
                        zip_data = crear_archivo_zip(archivos_xyz)
                        st.download_button(
                            label="📦 Descargar Archivos XYZ (ZIP)",
                            data=zip_data,
                            file_name="estereoisomeros_xyz.zip",
                            mime="application/zip",
                            key="download_zip_button"
                        )
                        with st.expander("👀 Vista Previa (Primer XYZ)"):
                            primer_archivo = list(archivos_xyz.values())[0]
                            st.code(primer_archivo)
        else:
            if smiles_input:
                 st.info("💡 Ingresa un SMILES con centros quirales especificados para iniciar la generación. Verifica el formato (@ o @@).")
        
        
    st.markdown("---")
    st.markdown(
        f"""
        <div style='text-align: center; color: #555555;'>
            <small>🧬 **Inchiral** | Desarrollado con Streamlit & RDKit | {BACKGROUND_COLOR} Theme</small>
        </div>
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()
