from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException, TimeoutException
from webdriver_manager.chrome import ChromeDriverManager
import time

# 设置Chrome驱动
options = Options()
# options.add_argument('--headless')  # 无头模式 - 先注释掉尝试非无头模式
options.add_argument('--disable-gpu')  # 禁用GPU加速
options.add_argument('--no-sandbox')  # 解决DevToolsActivePort文件不存在的报错
options.binary_location = "D:/Google/Chrome/Application/chrome.exe"  # 将此路径替换为你的Chrome安装路径

# 创建WebDriver实例
service = Service(ChromeDriverManager().install())
driver = webdriver.Chrome(service=service, options=options)

# 设置超时时间
driver.set_page_load_timeout(30)  # 30秒超时

# 请求初始网页
url = 'https://abrc.osu.edu/stocks/385307'  # 替换为目标网页的URL

try:
    driver.get(url)
    
    # 等待页面加载完成并获取标题
    WebDriverWait(driver, 10).until(EC.title_contains('ABRC'))

    # 打开文件以写入
    with open('output_links.txt', 'w', encoding='utf-8') as file:
        # 获取初始页面上的所有链接
        links = driver.find_elements(By.TAG_NAME, 'a')

        # 遍历所有链接
        for link in links:
            link_text = link.text
            link_url = link.get_attribute('href')
            if link_text.startswith('CS'):
                # 将链接文本和 URL 写入文件
                file.write(f"链接文本: {link_text}, URL: {link_url}\n")
                # 打印链接文本和 URL
                print(f"链接文本: {link_text}, URL: {link_url}")

        # 尝试点击从2到20的“加载更多”按钮
        for page_num in range(2, 21):  # 从2到20
            buttons_found = False
            try:
                # 尝试找到“加载更多”按钮
                load_more_button = driver.find_element(By.XPATH, f"//button[contains(text(), '{page_num}')]")
                if load_more_button:
                    load_more_button.click()
                    print(f"点击了'加载更多'按钮，第 {page_num} 页")
                    
                    # 等待页面加载新内容
                    WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.TAG_NAME, 'body')))
                    time.sleep(2)  # 等待新的内容加载

                    # 获取新加载内容的所有链接
                    new_links = driver.find_elements(By.TAG_NAME, 'a')
                    
                    # 检查是否有新的链接
                    new_links_found = False
                    for link in new_links:
                        link_text = link.text
                        link_url = link.get_attribute('href')
                        if link_text.startswith('CS'):
                            # 将链接文本和 URL 写入文件
                            file.write(f"链接文本: {link_text}, URL: {link_url}\n")
                            # 打印链接文本和 URL
                            print(f"链接文本: {link_text}, URL: {link_url}")
                            new_links_found = True
                    
                    if new_links_found:
                        buttons_found = True
                    else:
                        print(f"第 {page_num} 页没有新链接")
                    
            except NoSuchElementException:
                print(f"第 {page_num} 页的'加载更多'按钮不存在")
                
            # 仅在有新的按钮被点击时继续，否则退出循环
            if not buttons_found:
                print("没有找到更多的'加载更多'按钮，退出循环")
                break

except TimeoutException:
    print("页面加载超时")

except Exception as e:
    print(f"发生错误: {e}")

finally:
    driver.quit()

# 打开文件并读取内容
with open('output_links.txt', 'r', encoding='utf-8') as file:
    lines = file.readlines()

# 提取第四列，即 URL
fourth_column = [line.split()[-1].strip() for line in lines]

# 重新启动 WebDriver 实例
driver = webdriver.Chrome(service=service, options=options)

# 遍历所有链接 URL 并进行访问
for link_url in fourth_column:
    try:
        driver.get(link_url)
        # 等待页面加载完成并获取标题
        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.TAG_NAME, 'title')))
        
        # 打印访问的页面标题和URL
        print(f"访问的页面标题: {driver.title}, URL: {link_url}")
        
        # 在这里添加对每个链接页面的进一步处理代码，获取并打印页面文本
        body_text = driver.find_element(By.TAG_NAME, 'body').text

        # 将页面文本分割为行
        lines = body_text.split('\n')
        
        # 打开文件以写入模式
        with open('output_specie_3.txt', 'a', encoding='utf-8') as file:
            genome_id_line = None

            for i in range(len(lines)):

                # 查找以'1001 Genomes id:'开头的行
                if lines[i].strip().startswith('1001 Genomes id:'):
                    genome_id_line = lines[i].strip()
                    
                    # 从当前位置开始查找以'Species'开头的行及其下一行
                    for j in range(i + 1, len(lines)):
                        if lines[j].strip().startswith('Species'):
                            if j + 1 < len(lines):
                                species_line = lines[j + 1].strip()
                            break
                    
                    # 组合并写入文件
                    combined_line = f"{genome_id_line}, {species_line}\n"
                    file.write(combined_line)
                    break  # 找到后直接跳出循环

            # 如果没有找到 '1001 Genomes id:' 行，则查找 'Species' 行
            if not genome_id_line:
                for i in range(len(lines)):
                    if lines[i].strip().startswith('Species'):
                        if i + 1 < len(lines):
                            species_line = lines[i + 1].strip()
                            combined_line = f"Unknown 1001 Genomes id, {species_line}\n"
                            file.write(combined_line)
                            break

        
    except Exception as e:
        print(f"访问 {link_url} 时发生错误: {e}")

driver.quit()


