from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt

template = Path(r'C:\Users\26702\Desktop\9.23汇报.pptx')
out = Path(r'C:\Users\26702\Desktop\10.14汇报.pptx')
img_root = Path(r'E:\论文路径')

prs = Presentation(str(template))

def add_title_slide(title):
    slide = prs.slides.add_slide(prs.slide_layouts[0])
    slide.shapes.title.text = title

def add_image_slide(title, img):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    box = slide.shapes.add_textbox(Inches(0.4), Inches(0.15), Inches(12), Inches(0.4))
    box.text_frame.text = title
    box.text_frame.paragraphs[0].font.size = Pt(24)
    slide.shapes.add_picture(str(img), Inches(0.3), Inches(0.75), width=Inches(12.6))

add_title_slide('Carbon fixation pathway discovery and visualization\n10.14 update')

for module in ['EMP', 'ED', 'PKT', 'OPPP']:
    imgs = sorted((img_root / module).glob('*.png'))
    s = prs.slides.add_slide(prs.slide_layouts[1])
    s.shapes.title.text = module + ' module pathways'
    s.placeholders[1].text = f'{len(imgs)} candidate pathways rendered by the metabolic pathway visualization pipeline.'
    for img in imgs:
        add_image_slide(f'{module} | {img.stem}', img)

prs.save(str(out))
print(out)
