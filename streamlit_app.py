import streamlit as st

st.title("NeoApop-AI 🚀")
st.write("Salom! Bu sening birinchi AI + Drug Discovery dasturing!")

target = st.text_input("Kasallik yoki target protein nomini kiriting:")
if st.button("Molekulalar yaratish"):
    st.success(f"Demo: {target} uchun 1000 molekula yaratildi (test holat)")
import streamlit as st

st.title("NeoApop-AI 🚀")
st.write("Salom! Bu sening birinchi AI + Drug Discovery dasturing!")

target = st.text_input("Kasallik yoki target protein nomini kiriting:")
if st.button("Molekulalar yaratish"):
    st.success(f"Demo: {target} uchun 1000 molekula yaratildi (test holat)")
