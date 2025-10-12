# streamlit_app.py
import streamlit as st
import pandas as pd
import numpy as np
import hashlib
import requests
from io import BytesIO

# RDKit ixtiyoriy (agar o‘rnatilmagan bo‘lsa, app ishlashda davom etadi)
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors, Crippen, rdMolDescriptors
    RDKit_AVAILABLE = True
except ImportError:
    RDKit_AVAILABLE = False

# =============================
# 🔐 Streamlit Authenticator (new version)
# =============================
import streamlit_authenticator as stauth

# Foydalanuvchilar ro‘yxati
names = ["Admin"]
usernames = ["admin"]
passwords = ["demo123"]

# Cookie sozlamalari (yangi versiya sintaksisi)
cookie_settings = {
    "name": "neoapop_ai_cookie",
    "key": "neoapop_ai_sig",
    "expiry_days": 1
}

# Authenticator obyektini yaratish
authenticator = stauth.Authenticate(
    names=names,
    usernames=usernames,
    passwords=passwords,
    cookie_settings=cookie_settings
)

# Login oynasi
name, auth_status, username = authenticator.login("Login", "main")

# =============================
# 🌐 Agar login muvaffaqiyatli bo‘lsa
# =============================
if auth_status:
    st.sidebar.success(f"Xush kelibsiz, {name}!")
    authenticator.logout("Logout", "sidebar")

    st.title("NeoApop-AI 🚀 (Pro demo)")
    st.write("Salom! Bu kengaytirilgan AI + Drug Discovery platforma demo versiyasi.")

    # Tablar: Molekula / CSV / AI demo
    tab1, tab2, tab3 = st.tabs(["🧪 Molekula", "📁 CSV yuklash", "🧠 AI demo"])

    # ================
    # Tab 1 - Molekula tahlili
    # ================
    with tab1:
        st.header("🧬 Molekula tahlili (SMILES orqali)")
        smiles = st.text_input("SMILES kiriting (masalan: CC(=O)O yoki c1ccccc1):")

        if st.button("Tahlil qilish"):
            if not smiles:
                st.warning("Iltimos SMILES kiriting.")
            elif not RDKit_AVAILABLE:
                st.error("RDKit o‘rnatilmagan — descriptorlar hisoblanmaydi.")
            else:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    props = {
                        "MolWt": round(Descriptors.MolWt(mol), 3),
                        "LogP": round(Crippen.MolLogP(mol), 3),
                        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
                        "HBA": rdMolDescriptors.CalcNumHBA(mol),
                        "HBD": rdMolDescriptors.CalcNumHBD(mol),
                        "RotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol)
                    }
                    st.table(pd.DataFrame([props], index=[smiles]))
                    img = Draw.MolToImage(mol, size=(300, 250))
                    st.image(img, caption="2D strukturasi")
                else:
                    st.error("Noto‘g‘ri SMILES format.")
    
    # ================
    # Tab 2 - CSV yuklash
    # ================
    with tab2:
        st.header("📂 CSV yuklash (batch tahlil)")
        st.write("Faylda `smiles` nomli ustun bo‘lishi kerak.")
        file = st.file_uploader("CSV faylni tanlang", type=["csv"])

        if file is not None:
            df = pd.read_csv(file)
            st.write("Fayldagi ma'lumotlar:")
            st.dataframe(df.head())

            if "smiles" not in df.columns:
                st.error("CSV faylda `smiles` ustuni topilmadi.")
            elif st.button("Hisoblashni boshlash"):
                if not RDKit_AVAILABLE:
                    st.warning("RDKit mavjud emas, descriptorlar hisoblanmaydi.")
                else:
                    results = []
                    for smi in df["smiles"]:
                        mol = Chem.MolFromSmiles(smi)
                        if mol:
                            props = {
                                "SMILES": smi,
                                "MolWt": round(Descriptors.MolWt(mol), 3),
                                "LogP": round(Crippen.MolLogP(mol), 3),
                                "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3)
                            }
                            results.append(props)
                    res_df = pd.DataFrame(results)
                    st.success(f"{len(res_df)} ta molekula tahlil qilindi ✅")
                    st.dataframe(res_df.head())
                    csv = res_df.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 CSV natijani yuklash", csv, "results.csv")

    # ================
    # Tab 3 - AI demo
    # ================
    with tab3:
        st.header("🧠 AI faollik bashorati (demo)")
        st.write("Bu demo versiyada tasodifiy AI baholash ishlatiladi.")
        smi = st.text_input("SMILES kiriting (masalan: CCO):")
        if st.button("AI baholash"):
            if smi:
                score = round(np.random.uniform(40, 95), 2)
                st.metric("Taxminiy faollik", f"{score} %")
            else:
                st.warning("SMILES kiriting.")

# =============================
# 🚫 Login muvaffaqiyatsiz
# =============================
elif auth_status is False:
    st.error("Login yoki parol noto‘g‘ri ❌")
else:
    st.warning("Iltimos login va parolni kiriting.")
# streamlit_app.py
import streamlit as st
import pandas as pd
import numpy as np
import hashlib
import requests
from io import BytesIO

# RDKit ixtiyoriy (agar o‘rnatilmagan bo‘lsa, app ishlashda davom etadi)
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors, Crippen, rdMolDescriptors
    RDKit_AVAILABLE = True
except ImportError:
    RDKit_AVAILABLE = False

# =============================
# 🔐 Streamlit Authenticator (new version)
# =============================
import streamlit_authenticator as stauth

# Foydalanuvchilar ro‘yxati
names = ["Admin"]
usernames = ["admin"]
passwords = ["demo123"]

# Cookie sozlamalari (yangi versiya sintaksisi)
cookie_settings = {
    "name": "neoapop_ai_cookie",
    "key": "neoapop_ai_sig",
    "expiry_days": 1
}

# Authenticator obyektini yaratish
authenticator = stauth.Authenticate(
    names=names,
    usernames=usernames,
    passwords=passwords,
    cookie_settings=cookie_settings
)

# Login oynasi
name, auth_status, username = authenticator.login("Login", "main")

# =============================
# 🌐 Agar login muvaffaqiyatli bo‘lsa
# =============================
if auth_status:
    st.sidebar.success(f"Xush kelibsiz, {name}!")
    authenticator.logout("Logout", "sidebar")

    st.title("NeoApop-AI 🚀 (Pro demo)")
    st.write("Salom! Bu kengaytirilgan AI + Drug Discovery platforma demo versiyasi.")

    # Tablar: Molekula / CSV / AI demo
    tab1, tab2, tab3 = st.tabs(["🧪 Molekula", "📁 CSV yuklash", "🧠 AI demo"])

    # ================
    # Tab 1 - Molekula tahlili
    # ================
    with tab1:
        st.header("🧬 Molekula tahlili (SMILES orqali)")
        smiles = st.text_input("SMILES kiriting (masalan: CC(=O)O yoki c1ccccc1):")

        if st.button("Tahlil qilish"):
            if not smiles:
                st.warning("Iltimos SMILES kiriting.")
            elif not RDKit_AVAILABLE:
                st.error("RDKit o‘rnatilmagan — descriptorlar hisoblanmaydi.")
            else:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    props = {
                        "MolWt": round(Descriptors.MolWt(mol), 3),
                        "LogP": round(Crippen.MolLogP(mol), 3),
                        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
                        "HBA": rdMolDescriptors.CalcNumHBA(mol),
                        "HBD": rdMolDescriptors.CalcNumHBD(mol),
                        "RotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol)
                    }
                    st.table(pd.DataFrame([props], index=[smiles]))
                    img = Draw.MolToImage(mol, size=(300, 250))
                    st.image(img, caption="2D strukturasi")
                else:
                    st.error("Noto‘g‘ri SMILES format.")
    
    # ================
    # Tab 2 - CSV yuklash
    # ================
    with tab2:
        st.header("📂 CSV yuklash (batch tahlil)")
        st.write("Faylda `smiles` nomli ustun bo‘lishi kerak.")
        file = st.file_uploader("CSV faylni tanlang", type=["csv"])

        if file is not None:
            df = pd.read_csv(file)
            st.write("Fayldagi ma'lumotlar:")
            st.dataframe(df.head())

            if "smiles" not in df.columns:
                st.error("CSV faylda `smiles` ustuni topilmadi.")
            elif st.button("Hisoblashni boshlash"):
                if not RDKit_AVAILABLE:
                    st.warning("RDKit mavjud emas, descriptorlar hisoblanmaydi.")
                else:
                    results = []
                    for smi in df["smiles"]:
                        mol = Chem.MolFromSmiles(smi)
                        if mol:
                            props = {
                                "SMILES": smi,
                                "MolWt": round(Descriptors.MolWt(mol), 3),
                                "LogP": round(Crippen.MolLogP(mol), 3),
                                "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3)
                            }
                            results.append(props)
                    res_df = pd.DataFrame(results)
                    st.success(f"{len(res_df)} ta molekula tahlil qilindi ✅")
                    st.dataframe(res_df.head())
                    csv = res_df.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 CSV natijani yuklash", csv, "results.csv")

    # ================
    # Tab 3 - AI demo
    # ================
    with tab3:
        st.header("🧠 AI faollik bashorati (demo)")
        st.write("Bu demo versiyada tasodifiy AI baholash ishlatiladi.")
        smi = st.text_input("SMILES kiriting (masalan: CCO):")
        if st.button("AI baholash"):
            if smi:
                score = round(np.random.uniform(40, 95), 2)
                st.metric("Taxminiy faollik", f"{score} %")
            else:
                st.warning("SMILES kiriting.")

# =============================
# 🚫 Login muvaffaqiyatsiz
# =============================
elif auth_status is False:
    st.error("Login yoki parol noto‘g‘ri ❌")
else:
    st.warning("Iltimos login va parolni kiriting.")
# streamlit_app.py
import streamlit as st
import pandas as pd
import numpy as np
import hashlib
import requests
from io import BytesIO

# RDKit ixtiyoriy (agar o‘rnatilmagan bo‘lsa, app ishlashda davom etadi)
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors, Crippen, rdMolDescriptors
    RDKit_AVAILABLE = True
except ImportError:
    RDKit_AVAILABLE = False

# =============================
# 🔐 Streamlit Authenticator (new version)
# =============================
import streamlit_authenticator as stauth

# Foydalanuvchilar ro‘yxati
names = ["Admin"]
usernames = ["admin"]
passwords = ["demo123"]

# Cookie sozlamalari (yangi versiya sintaksisi)
cookie_settings = {
    "name": "neoapop_ai_cookie",
    "key": "neoapop_ai_sig",
    "expiry_days": 1
}

# Authenticator obyektini yaratish
authenticator = stauth.Authenticate(
    names=names,
    usernames=usernames,
    passwords=passwords,
    cookie_settings=cookie_settings
)

# Login oynasi
name, auth_status, username = authenticator.login("Login", "main")

# =============================
# 🌐 Agar login muvaffaqiyatli bo‘lsa
# =============================
if auth_status:
    st.sidebar.success(f"Xush kelibsiz, {name}!")
    authenticator.logout("Logout", "sidebar")

    st.title("NeoApop-AI 🚀 (Pro demo)")
    st.write("Salom! Bu kengaytirilgan AI + Drug Discovery platforma demo versiyasi.")

    # Tablar: Molekula / CSV / AI demo
    tab1, tab2, tab3 = st.tabs(["🧪 Molekula", "📁 CSV yuklash", "🧠 AI demo"])

    # ================
    # Tab 1 - Molekula tahlili
    # ================
    with tab1:
        st.header("🧬 Molekula tahlili (SMILES orqali)")
        smiles = st.text_input("SMILES kiriting (masalan: CC(=O)O yoki c1ccccc1):")

        if st.button("Tahlil qilish"):
            if not smiles:
                st.warning("Iltimos SMILES kiriting.")
            elif not RDKit_AVAILABLE:
                st.error("RDKit o‘rnatilmagan — descriptorlar hisoblanmaydi.")
            else:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    props = {
                        "MolWt": round(Descriptors.MolWt(mol), 3),
                        "LogP": round(Crippen.MolLogP(mol), 3),
                        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
                        "HBA": rdMolDescriptors.CalcNumHBA(mol),
                        "HBD": rdMolDescriptors.CalcNumHBD(mol),
                        "RotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol)
                    }
                    st.table(pd.DataFrame([props], index=[smiles]))
                    img = Draw.MolToImage(mol, size=(300, 250))
                    st.image(img, caption="2D strukturasi")
                else:
                    st.error("Noto‘g‘ri SMILES format.")
    
    # ================
    # Tab 2 - CSV yuklash
    # ================
    with tab2:
        st.header("📂 CSV yuklash (batch tahlil)")
        st.write("Faylda `smiles` nomli ustun bo‘lishi kerak.")
        file = st.file_uploader("CSV faylni tanlang", type=["csv"])

        if file is not None:
            df = pd.read_csv(file)
            st.write("Fayldagi ma'lumotlar:")
            st.dataframe(df.head())

            if "smiles" not in df.columns:
                st.error("CSV faylda `smiles` ustuni topilmadi.")
            elif st.button("Hisoblashni boshlash"):
                if not RDKit_AVAILABLE:
                    st.warning("RDKit mavjud emas, descriptorlar hisoblanmaydi.")
                else:
                    results = []
                    for smi in df["smiles"]:
                        mol = Chem.MolFromSmiles(smi)
                        if mol:
                            props = {
                                "SMILES": smi,
                                "MolWt": round(Descriptors.MolWt(mol), 3),
                                "LogP": round(Crippen.MolLogP(mol), 3),
                                "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3)
                            }
                            results.append(props)
                    res_df = pd.DataFrame(results)
                    st.success(f"{len(res_df)} ta molekula tahlil qilindi ✅")
                    st.dataframe(res_df.head())
                    csv = res_df.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 CSV natijani yuklash", csv, "results.csv")

    # ================
    # Tab 3 - AI demo
    # ================
    with tab3:
        st.header("🧠 AI faollik bashorati (demo)")
        st.write("Bu demo versiyada tasodifiy AI baholash ishlatiladi.")
        smi = st.text_input("SMILES kiriting (masalan: CCO):")
        if st.button("AI baholash"):
            if smi:
                score = round(np.random.uniform(40, 95), 2)
                st.metric("Taxminiy faollik", f"{score} %")
            else:
                st.warning("SMILES kiriting.")

# =============================
# 🚫 Login muvaffaqiyatsiz
# =============================
elif auth_status is False:
    st.error("Login yoki parol noto‘g‘ri ❌")
else:
    st.warning("Iltimos login va parolni kiriting.")


