import cv2
import os
import re

def create_video_from_images(image_folder, video_name, fps):
    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape
    pattern = re.compile(r'\d+')

    # 对文件名中的数字进行提取和排序
    images.sort(key=lambda x: int(re.search(pattern, x).group()))

    video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), fps, (width,height))

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    cv2.destroyAllWindows()
    video.release()

# 使用示例
image_folder = 'KB/fac=0.2_dt=0.01_nu=0.0002_pert=0.01_vi=1' # 图片所在文件夹的路径
video_name = 'KB_output_video.mp4' # 视频输出文件名
fps = 12 # 视频帧率

create_video_from_images(image_folder, video_name, fps)