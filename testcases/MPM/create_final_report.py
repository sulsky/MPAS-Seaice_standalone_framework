from pypdf import PdfReader, PdfWriter, PageObject, PdfWriter
import img2pdf
import os

#-------------------------------------------------------------------------------

def stitch_pdfs_side_by_side(pdf1_path,
                             pdf2_path,
                             output_path):

    reader1 = PdfReader(pdf1_path)
    reader2 = PdfReader(pdf2_path)
    writer = PdfWriter()

    # Grab the first page from each PDF
    page1 = reader1.pages[0]
    page2 = reader2.pages[0]

    # Extract dimensions
    width1 = page1.mediabox.width
    height1 = page1.mediabox.height
    width2 = page2.mediabox.width
    height2 = page2.mediabox.height

    # Calculate dimensions for the wide canvas
    total_width = width1 + width2
    max_height = max(height1, height2)

    # Create the modern blank layout page
    new_page = PageObject.create_blank_page(width=total_width, height=max_height)

    # Merge pages onto the canvas
    new_page.merge_page(page1)
    new_page.merge_translated_page(page2, tx=width1, ty=0)

    # Write out the final canvas
    writer.add_page(new_page)
    with open(output_path, "wb") as f:
        writer.write(f)

#-------------------------------------------------------------------------------

def create_final_report():

    """
    This function combines testcase plots from this checkout of the testcases with
    a baseline checkout of the testcases
    """

    MPAS_SEAICE_TESTCASE_BASELINES = os.environ.get('MPAS_SEAICE_TESTCASE_BASELINES')
    if (MPAS_SEAICE_TESTCASE_BASELINES is None):
        raise Exception("MPAS_SEAICE_TESTCASE_BASELINES must be set")

    testLocation = "../../"

    plots = ["/testcases/MPM/polynya/cells_scatter.png",
             "/testcases/MPM/polynya/particles_scatter.png",
             "/testcases/MPM/spherical_operators/interpolation/interpolation_scaling.png",
             "/testcases/MPM/spherical_operators/reconstruction/reconstruction_map.png",
             "/testcases/MPM/spherical_operators/reconstruction/reconstruction_scaling.png"]

    writer = PdfWriter()

    for i, plot in enumerate(plots):

        filename1 = MPAS_SEAICE_TESTCASE_BASELINES+plot
        filename2 = testLocation+plot

        with open("tmp1.pdf", "wb") as f:
            f.write(img2pdf.convert(filename1))
        with open("tmp2.pdf", "wb") as f:
            f.write(img2pdf.convert(filename2))
        stitch_pdfs_side_by_side("tmp1.pdf", "tmp2.pdf", "tmp3.pdf")

        writer.append("tmp3.pdf")

    with open("final_report.pdf", "wb") as f:
        writer.write(f)

    writer.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_final_report()
