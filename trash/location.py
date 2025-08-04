from geopy.geocoders import GoogleV3

# 创建一个 GoogleV3 geolocator 对象，设置你的 API 密钥
geolocator = GoogleV3(api_key='YOUR_API_KEY')

# 查询国家名称
location = geolocator.geocode("Armenia")

# 打印经纬度
if location:
    print(f"Latitude: {location.latitude}, Longitude: {location.longitude}")
else:
    print("Location not found")
