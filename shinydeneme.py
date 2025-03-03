from shiny import App, ui, render
import zipfile
import os 
import shutil #dosyayı kopyalama-taşıma
import scanpy as sc
import matplotlib.pyplot as plt

UPLOAD_FOLDER = "uploads" 
os.makedirs(UPLOAD_FOLDER, exist_ok=True) #klasör yoksa otomatik olarak oluşturmak için

def process_h5ad(h5ad_path):
    """H5AD dosyasını işler ve UMAP görselleştirmesi yapar."""
    adata = sc.read_h5ad(h5ad_path) 
    sc.pp.pca(adata) 
    sc.pp.neighbors(adata) 
    sc.tl.umap(adata) #umap hesaplar
    sc.pl.umap(adata, #umap görselleştirir
               color=adata.obs.columns[0],
               legend_loc="right", #görseldeki renklerin görselin sağında çıkması için
               legend_fontsize=10,
               show=False,
    )  
    plt.savefig("umap.png")

# shıny nin arayüzünü oluşturur.
app_ui = ui.page_fluid(
    ui.panel_title("ZIP Dosyası İşleyici ve UMAP Görselleştirici"),
    ui.input_file("zip_file", "Lütfen bir ZIP dosyası yükleyin", accept=[".zip"], multiple=False),
    ui.output_text("status"),
    ui.output_image("umap_plot")
)

def server(input, output, session):
    @output
    @render.text
    def status(): # zip dosyasını yükler, açar, işler.
        file = input.zip_file() 
        if not file:
            return "Henüz bir dosya yüklenmedi."

        # Dosya bilgilerini upload kalsörüne kopyalar (shuthil fonk.).
        file_info = file[0]  # İlk dosyayı al
        file_path = os.path.join(UPLOAD_FOLDER, file_info["name"])
        shutil.copy(file_info["datapath"], file_path)

        # ZIP dosyasını aç
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(UPLOAD_FOLDER)

        # process_h5ad fonk. .h5ad dosyasını bulur, işler.
        h5ad_file = [f for f in os.listdir(UPLOAD_FOLDER) if f.endswith(".h5ad")][0]
        process_h5ad(os.path.join(UPLOAD_FOLDER, h5ad_file))
        return "ZIP dosyası başarıyla yüklendi ve UMAP görselleştirildi!"

    @output
    @render.image
    def umap_plot():
        if not input.zip_file():
            return None
        return {"src": "umap.png"}  # görselin pathini shiny e iletir.

app = App(app_ui, server)
