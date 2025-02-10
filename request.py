from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

# 设置Chrome驱动
options = Options()
options.add_argument('--headless')  # 如果需要无头模式
service = Service(ChromeDriverManager().install())
driver = webdriver.Chrome(service=service, options=options)

# 请求网页
url = 'https://abrc.osu.edu/stocks/322498'  # 请替换为实际URL
driver.get(url)

try:
    # 等待特定元素加载完成
    WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.LINK_TEXT, "CS76347")))

    # 查找特定链接
    link_element = driver.find_element(By.LINK_TEXT, "CS76347")
    specific_link = link_element.get_attribute('href')

    print(f"The link for 'CS76347' is: {specific_link}")
finally:
    driver.quit()
