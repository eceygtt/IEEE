

UPLOAD_FOLDER = "uploads"
STATIC_FOLDER = "static"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(STATIC_FOLDER, exist_ok=True)

def unzip_file(zip_path, extract_to):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)
    h5ad_files = [f for f in os.listdir(extract_to) if f.endswith(".h5ad")]
    if not h5ad_files:
        raise FileNotFoundError("ZIP dosyasında .h5ad dosyası bulunamadı.")
    return os.path.join(extract_to, h5ad_files[0])

# Function to process the .h5ad file and generate UMAP plot
def process_h5ad(h5ad_path, output_folder):
    adata = sc.read_h5ad(h5ad_path)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    
    # Plot UMAP
    plt.figure(figsize=(8, 6))
    sc.pl.umap(adata, color="louvain", show=False)
    plot_path = os.path.join(output_folder, "umap_plot.png")
    plt.savefig(plot_path)
    plt.close()
    return plot_path

# Define the UI
app_ui = ui.page_fluid(
    ui.panel_title("Tek Hücre RNA-seq UMAP Görselleştirici"),
    ui.input_file("zip_file", "Lütfen ZIP dosyanızı buraya sürükleyin veya seçin", accept=[".zip"]),
    ui.output_text("file_status"),
    ui.output_image("umap_plot")
)

# Define the server logic
def server(input, output, session):
    @output
    @render.text
    def file_status():
        file = input.zip_file()
        if file is None:
            return "Henüz bir dosya yüklenmedi."
        
        # Save the uploaded file
        file_path = os.path.join(UPLOAD_FOLDER, file.name)
        shutil.move(file["datapath"], file_path)
        
        try:
            # Unzip and process the file
            h5ad_path = unzip_file(file_path, UPLOAD_FOLDER)
            plot_path = process_h5ad(h5ad_path, STATIC_FOLDER)
            return "Dosya başarıyla işlendi ve UMAP görselleştirildi!"
        except Exception as e:
            return f"Hata: {str(e)}"
    
    @output
    @render.image
    def umap_plot():
        file = input.zip_file()
        if file is None:
            return None
        plot_path = os.path.join(STATIC_FOLDER, "umap_plot.png")
        if os.path.exists(plot_path):
            return {"src": plot_path, "width": "600px", "height": "450px"}
        return None

# Create the Shiny app
app = App(app_ui, server)
