# streamlit_app.py
"""
NeoApop-AI (Pro demo)
Single-file Streamlit app with:
- streamlit-authenticator login (demo hashed passwords)
- SMILES input (2D image shown via RDKit if available)
- Molecule descriptors (MW, LogP, TPSA, HBD/HBA, rotatable bonds)
- CSV upload (SMILES column)
- Snowflake skeleton example (commented)
- Simple AI prediction stub (random for demo)
- Graceful handling when RDKit or py3Dmol unavailable
"""

from typing import Optional
import streamlit as st
import pandas as pd
import numpy as np
import hashlib
import requests
from io import BytesIO

# Try to import RDKit and 3D renderer; if not available, disable those features gracefully
RDKit_AVAILABLE = True
TRY_3DMOL = True
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem, Descriptors, Crippen, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
except Exception as e:
    RDKit_AVAILABLE = False

try:
    import py3Dmol
    from stmol import showmol
except Exception:
    TRY_3DMOL = False

# Authenticator (streamlit-authenticator)
import streamlit_authenticator as stauth

# ------------------------------
# Helper utilities
# ------------------------------
def hash_password(password: str) -> str:
    return hashlib.sha256(password.encode()).hexdigest()

def descriptors_from_smiles(smiles: str) -> Optional[dict]:
    if not RDKit_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    props = {
        "SMILES": smiles,
        "MolWt": round(Descriptors.MolWt(mol), 3),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "HBA": int(rdMolDescriptors.CalcNumHBA(mol)),
        "HBD": int(rdMolDescriptors.CalcNumHBD(mol)),
        "RotatableBonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol))
    }
    return props

def mol_to_png(smiles: str, size=(300, 220)) -> Optional[BytesIO]:
    if not RDKit_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=size)
    buf = BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return buf

def try_embed_3d(smiles: str, height=350, width=350):
    if not (RDKit_AVAILABLE and TRY_3DMOL):
        st.info("3D viewer mavjud emas (serverda py3Dmol yoki RDKit o'rnatilmagan).")
        return
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    mb = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(mb, 'mol')
    viewer.setStyle({ 'stick':{} })
    viewer.zoomTo()
    showmol(viewer, height=height)

# ------------------------------
# Auth setup (demo)
# ------------------------------
# NOTE: In production use hashed passwords and env vars.
names = ["Admin"]
usernames = ["admin"]
# Example hashed passwords (demo only). Use stauth.Hasher([...]).generate() locally to create.
# Here we use plain text list but stauth will hash them internally. For demo it's OK.
passwords = ["demo123"]  

authenticator = stauth.Authenticate(
    names, usernames, passwords,
    "neoapop_ai_cookie", "neoapop_ai_sig", cookie_expiry_days=1
)

name, auth_status, username = authenticator.login("Login", "main")

if auth_status:
    st.sidebar.success(f"Xush kelibsiz, {name}!")

    # Sidebar settings
    st.sidebar.header("Settings")
    bg = st.sidebar.selectbox("Background", ["white","black","transparent"])
    language = st.sidebar.selectbox("Language", ["uz","en","ru"])
    show_3d_opt = st.sidebar.checkbox("Enable 3D viewer (if available)", value=True)

    # Main UI
    st.set_page_config(page_title="NeoApop-AI", layout="centered")
    st.title("NeoApop-AI 🚀 (Pro demo)")
    st.write("Salom! Bu NeoApop-AI platforma ning kengaytirilgan demo versiyasi.")

    # Tabs for organization
    tab1, tab2, tab3 = st.tabs(["Molecule", "Batch (CSV)", "Snowflake / AI"])

    # ---------------------
    # Tab 1: Single molecule
    # ---------------------
    with tab1:
        st.header("Single molecule analysis")
        col1, col2 = st.columns([2,1])
        with col1:
            smiles = st.text_input("SMILES kiritish (mas. CC(=O)O, c1ccccc1):")
            if st.button("Analyze molecule"):
                if not smiles:
                    st.warning("Iltimos SMILES kiriting.")
                else:
                    if RDKit_AVAILABLE:
                        props = descriptors_from_smiles(smiles)
                        if props is None:
                            st.error("Noto'g'ri SMILES yoki RDKit muammosi.")
                        else:
                            st.subheader("Molecule descriptors")
                            st.table(pd.DataFrame([props]).set_index("SMILES").T)
                            # 2D image
                            img_buf = mol_to_png(smiles, size=(400,300))
                            if img_buf:
                                st.image(img_buf, caption="2D structure (RDKit)")
                            else:
                                st.info("2D rasm hosil bo'lmadi (RDKit yo'q).")
                            # 3D viewer
                            if show_3d_opt and TRY_3DMOL and RDKit_AVAILABLE:
                                try:
                                    try_embed_3d(smiles, height=350, width=350)
                                except Exception as e:
                                    st.warning("3D vizualizatsiya ishga tushmadi.")
                            else:
                                if show_3d_opt:
                                    st.info("3D viewer mavjud emas (py3Dmol yoki RDKit yo'q).")
                    else:
                        st.warning("RDKit serverda mavjud emas — 2D/3D va descriptors o‘rnatilmagan.")
                        st.info("Lekin siz CSV upload orqali massaviy scout qilishingiz mumkin.")

        with col2:
            st.info("Molekulani tahlil qilish (descriptor + 2D/3D).")
            st.write("Pro tip: to'g'ri natija uchun valid SMILES kiriting.")
            if st.button("Random demo molecule"):
                demo = "CCO"
                st.write("Demo SMILES:", demo)
                props = descriptors_from_smiles(demo) if RDKit_AVAILABLE else None
                if props:
                    st.table(pd.DataFrame([props]).set_index("SMILES").T)
                    st.image(mol_to_png(demo), caption="Ethanol (2D)")
                else:
                    st.write("RDKit mavjud emas, demo faqat xabar sifatida.")

    # ---------------------
    # Tab 2: Batch upload
    # ---------------------
    with tab2:
        st.header("Batch processing (CSV)")
        st.write("CSV fayl yuklang. Faylda `smiles` nomli ustun bo'lishi kerak.")
        uploaded = st.file_uploader("CSV fayl yuklash (.csv)", type=["csv"])
        if uploaded:
            df = pd.read_csv(uploaded)
            st.write("Fayl namunasi:")
            st.dataframe(df.head())
            if "smiles" not in df.columns:
                st.error("CSV ichida 'smiles' nomli ustun topilmadi. Ustun nomini tekshiring.")
            else:
                if st.button("Compute descriptors for all"):
                    if not RDKit_AVAILABLE:
                        st.error("RDKit o'rnatilmagan — descriptors hisoblanmaydi.")
                    else:
                        results = []
                        for s in df['smiles'].astype(str):
                            p = descriptors_from_smiles(s)
                            if p:
                                results.append(p)
                        res_df = pd.DataFrame(results)
                        st.success(f"Descriptors hisoblandi: {len(res_df)} ta.")
                        st.dataframe(res_df.head())
                        csv = res_df.to_csv(index=False).encode('utf-8')
                        st.download_button("Download results CSV", csv, file_name="neoapop_descriptors.csv")

    # ---------------------
    # Tab 3: Snowflake skeleton + AI stub
    # ---------------------
    with tab3:
        st.header("Snowflake integration & AI prediction (skeleton)")
        st.write("Bu bo'limda haqiqiy Snowflake ulanishi va AI model natijalari ko'rsatilishi mumkin.")
        st.markdown("""
        **Snowflake**: agar ulanish kerak bo'lsa, `snowflake_config.py` ichida maxfiy ma'lumotlarni saqlang va quyidagi kodni yoqing:
        ```py
        from snowflake.snowpark import Session
        session = Session.builder.configs(SNOWFLAKE_CONFIG).create()
        df = session.sql("SELECT * FROM MY_TABLE").to_pandas()
        ```
        """)
        st.info("AI prediction bu demo uchun random ball beradi (production: o'zingizning modelni joylang).")
        demo_smiles = st.text_input("AI demo uchun SMILES (optional):")
        if st.button("Predict activity (demo)"):
            if not demo_smiles:
                st.warning("Iltimos SMILES kiriting.")
            else:
                # Demo: fake prediction using simple rule or random
                if RDKit_AVAILABLE:
                    mol = Chem.MolFromSmiles(demo_smiles)
                    if mol is None:
                        st.error("Noto'g'ri SMILES")
                    else:
                        # simple heuristic: higher MW -> lower "permeability" (fake)
                        mw = Descriptors.MolWt(mol)
                        score = max(0, 100 - mw/2 + np.random.normal(0,5))
                        st.metric("Predicted activity (demo score)", f"{score:.1f}")
                else:
                    st.metric("Predicted activity (demo score)", f"{np.random.uniform(20,80):.1f}")

    # Logout
    if st.sidebar.button("Logout"):
        authenticator.logout("main")
        st.experimental_rerun()

elif auth_status is False:
    st.error("Login yoki parol xato ❌")
else:
    st.warning("Iltimos login va parolni kiriting.")
