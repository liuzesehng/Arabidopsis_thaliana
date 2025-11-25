# conda create -n r_geo -c conda-forge r-base=4.3 r-terra r-geodata
library(geodata)
library(terra)

# 1. 下载气候数据 (如果 data 目录下已有文件，不会重复下载)
# 注意: res=10 是最低分辨率(约18.5km)，如果需要更高精度，可以改为 5, 2.5 或 0.5
clim_data <- worldclim_global(var = "bio", res = 0.5, path = "./data")

# 定义输入文件列表
input_files <- c("Ala.ab.location.txt", "Ala.normal.location.txt")

# 定义英文列名列表 (按顺序 1-19)
bio_names <- c(
  "Annual Mean Temp", "Mean Diurnal Range", "Isothermality", "Temp Seasonality",
  "Max Temp Warmest Month", "Min Temp Coldest Month", "Temp Annual Range",
  "Mean Temp Wettest Q", "Mean Temp Driest Q", "Mean Temp Warmest Q", "Mean Temp Coldest Q",
  "Annual Precip", "Precip Wettest Month", "Precip Driest Month", "Precip Seasonality",
  "Precip Wettest Q", "Precip Driest Q", "Precip Warmest Q", "Precip Coldest Q"
)

for (file_path in input_files) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  
  message(paste("Processing:", file_path))
  
  # 读取数据，假设是制表符分隔
  points <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  
  # 确保 lat 和 long 列存在且为数值型
  if (!all(c("lat", "long") %in% colnames(points))) {
    warning(paste("Skipping", file_path, "- missing lat/long columns"))
    next
  }
  
  # 转换坐标列为数值，处理非数值数据
  points$lat <- as.numeric(points$lat)
  points$long <- as.numeric(points$long)
  
  # 提取有效坐标行进行查询
  valid_rows <- !is.na(points$lat) & !is.na(points$long)
  
  if (sum(valid_rows) > 0) {
    coords <- points[valid_rows, c("long", "lat")]
    clim_values <- terra::extract(clim_data, coords, ID = FALSE)
    colnames(clim_values) <- bio_names
    
    # 将结果合并回原数据框
    # 初始化气候列为 NA
    for (col in bio_names) {
      points[[col]] <- NA
    }
    
    points[valid_rows, bio_names] <- clim_values
  } else {
    warning(paste("No valid coordinates found in", file_path))
  }
  
  # 生成输出文件名
  output_file <- sub(".txt$", ".climate.txt", file_path)
  
  # 保存结果
  write.table(points, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste("Saved to:", output_file))
}
