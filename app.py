import streamlit as st
import itertools
import os
import zipfile
import tempfile
import io
import sys
import warnings

# Evitar warnings de RDKit
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

# Importar py3Dmol para visualización 3D
import py3Dmol

def mostrar_molecula_xyz(xyz_string):
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(xyz_string, "xyz")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer

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
                posiciones.append((i, True))  # ya es @@
                i += 2
            else:
                posiciones.append((i, False))  # es @ simple
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

def main():
    st.set_page_config(
        page_title="Inchiral - Generador de Estereoisómeros",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.title("🧬 Generador de Estereoisómeros")
    st.markdown("**Genera todos los estereoisómeros posibles y convierte a formato XYZ con visualización 3D**")
    
    with st.sidebar:
        try:
            st.image("imagenes1/logo.png", width=200)
        except:
            st.markdown("**🧬 Inchiral**")
        
        st.markdown("---")
        st.title("ℹ️ Información")
        st.markdown("""
        **Instrucciones:**
        1. Escribe un código SMILES, ya sea con quiralidad definida o no.
        2. La aplicación identificará de forma automática si la molécula presenta quiralidad.
        3. Cuando el SMILES incluye centros quirales indicados (@ o @@), se generarán todos los estereoisómeros posibles.
        4. Se admite un máximo de 3 centros quirales para evitar un número excesivo de isómeros.
        5. Opcionalmente convierte a formato XYZ y visualiza la molécula en 3D.
        
        **Ejemplos de SMILES:**
        - Molécula simple (sin centros quirales): C=CC
        - Molécula con un centro quiral: CC(O)C
        - Con estereoquímica definida: C[C@H](N)O
        - Aminoácido alternativo (alanina): N[C@H](C)C(=O)O
        """)
    
    st.subheader("📝 Entrada de Datos")
    smiles_input = st.text_input(
        "👉 Ingresa el código SMILES:",
        placeholder="Ejemplo: C[C@H](O)[C@@H](N)C"
    )
    
    if smiles_input:
        st.subheader("🔍 Análisis de Quiralidad")
        
        es_quiral, mensaje_quiralidad, centros_detectados = detectar_quiralidad(smiles_input)
        centros_especificados, posiciones_at = analizar_centros_existentes(smiles_input)

        col1, col2 = st.columns(2)
        
        with col1:
            st.info(f"**🔎 Análisis con RDKit:**")
            if RDKIT_AVAILABLE:
                if es_quiral:
                    st.success(f"✅ {mensaje_quiralidad}")
                    if centros_detectados:
                        st.write("**Centros detectados:**")
                        for i, (idx, chirality) in enumerate(centros_detectados):
                            tipo_quiralidad = str(chirality) if chirality else "Sin asignar"
                            st.write(f"• Átomo {idx}: {tipo_quiralidad}")
                else:
                    if "inválido" in mensaje_quiralidad:
                        st.error(f"❌ {mensaje_quiralidad}")
                    else:
                        st.warning(f"⚠️ {mensaje_quiralidad}")
            else:
                st.warning("⚠️ RDKit no disponible para análisis")
        
        with col2:
            st.info(f"**📋 Centros especificados en SMILES:**")
            if centros_especificados > 0:
                st.success(f"✅ {centros_especificados} centros con @ o @@ especificados")
                for pos in posiciones_at:
                    st.write(f"• Posición {pos}")
            else:
                st.warning("⚠️ No hay centros especificados con @ o @@")
        
        if RDKIT_AVAILABLE and es_quiral and centros_especificados == 0:
            st.info("""
            💡 Tu molécula es quiral pero no tiene centros especificados con @ o @@.
            Ejemplo: `CC(O)C(N)C` → `C[C@H](O)[C@@H](N)C`
            """)
        
        # Inicializar
        isomeros, n_centros = [], 0
        if centros_especificados > 0:
            with st.spinner("🔄 Generando estereoisómeros..."):
                isomeros, n_centros = generar_estereoisomeros(smiles_input)
        
        tab1, tab2, tab3 = st.tabs(["📋 Lista Completa", "💾 Descargar SMI", "🧪 Convertir a XYZ"])
        
        if isomeros:
            with tab1:
                col1, col2 = st.columns(2)
                for i, isomero in enumerate(isomeros):
                    if i % 2 == 0:
                        col1.code(f"{i+1}. {isomero}")
                    else:
                        col2.code(f"{i+1}. {isomero}")
            
            with tab2:
                smi_content = "\n".join(isomeros)
                st.download_button(
                    label="📥 Descargar archivo.smi",
                    data=smi_content,
                    file_name="estereoisomeros.smi",
                    mime="text/plain"
                )
                with st.expander("👀 Vista previa del archivo SMI"):
                    st.text(smi_content)
            
            with tab3:
                if st.button("🚀 Convertir todos a XYZ", type="primary"):
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    archivos_xyz = {}
                    mensajes_log = []
                    
                    for i, smiles in enumerate(isomeros):
                        try:
                            progress = (i + 1) / len(isomeros)
                            progress_bar.progress(progress)
                            status_text.text(f"Procesando molécula {i+1}/{len(isomeros)}: {smiles}")
                            
                            xyz_content, mensaje = smiles_to_xyz(smiles, i+1)
                            mensajes_log.append(mensaje)
                            
                            if xyz_content:
                                archivos_xyz[f"mol_{i+1}.xyz"] = xyz_content
                        except Exception as e:
                            mensajes_log.append(f"❌ Error procesando molécula {i+1}: {str(e)}")
                    
                    progress_bar.progress(1.0)
                    status_text.text("✅ Proceso completado!")
                    
                    with st.expander("📋 Log de procesamiento"):
                        for mensaje in mensajes_log:
                            if "❌" in mensaje or "⚠️" in mensaje:
                                st.error(mensaje)
                            else:
                                st.success(mensaje)
                    
                    if archivos_xyz:
                        zip_data = crear_archivo_zip(archivos_xyz)
                        st.download_button(
                            label="📦 Descargar archivos XYZ (ZIP)",
                            data=zip_data,
                            file_name="estereoisomeros_xyz.zip",
                            mime="application/zip"
                        )

                        # Selección para vista 3D
                        opciones = list(archivos_xyz.keys())
                        seleccion = st.selectbox("🧩 Selecciona un isómero para visualizar en 3D", opciones)
                        if seleccion:
                            st.subheader(f"Vista 3D de {seleccion}")
                            xyz_data = archivos_xyz[seleccion]
                            viewer = mostrar_molecula_xyz(xyz_data)
                            html = viewer._make_html()
                            st.components.v1.html(html, height=450)
        else:
            st.info("💡 Ingresa un SMILES con centros quirales especificados (@ o @@) para generar estereoisómeros")
    
    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center'>
            <small>🧬 <strong>Inchiral</strong> - Universidad Científica del Sur<br>
            Generador de Estereoisómeros | Desarrollado con Streamlit, RDKit y py3Dmol</small>
        </div>
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()
