"""Reusable formal landscape-PDF primitives for retained DynSysLab reports."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle


@dataclass(frozen=True)
class FormalReportTheme:
    navy: str = "#17212b"
    blue: str = "#2a6f9b"
    pale: str = "#f6f8fa"
    border: str = "#9aa8b3"
    muted: str = "#5b6670"
    light_text: str = "#e7edf2"


class FormalLandscapePDF:
    """Small repeatable presentation layer used by report exporters."""

    def __init__(self, output: str | Path, *, title: str, footer_label: str,
                 theme: FormalReportTheme | None = None):
        self.output = Path(output).resolve()
        self.output.parent.mkdir(parents=True, exist_ok=True)
        self.theme = theme or FormalReportTheme()
        self.width, self.height = landscape(letter)
        self.canvas = canvas.Canvas(
            str(self.output), pagesize=(self.width, self.height), pageCompression=1
        )
        self.canvas.setTitle(title)
        self.canvas.setAuthor("DynSysLab")
        self.footer_label = footer_label
        self.page_number = 1
        self.styles = getSampleStyleSheet()

    def footer(self) -> None:
        pdf = self.canvas
        pdf.setStrokeColor(colors.HexColor("#c5ced6"))
        pdf.line(30, 25, self.width - 30, 25)
        pdf.setFillColor(colors.HexColor(self.theme.muted))
        pdf.setFont("Helvetica", 7.5)
        pdf.drawString(34, 13, self.footer_label)
        pdf.drawRightString(self.width - 34, 13, f"Page {self.page_number}")

    def finish_page(self) -> None:
        self.footer()
        self.canvas.showPage()
        self.page_number += 1

    def title(self, text: str, *, size: float = 18.0) -> None:
        self.canvas.setFillColor(colors.HexColor(self.theme.navy))
        self.canvas.setFont("Helvetica-Bold", size)
        self.canvas.drawString(38, self.height - 40, text)

    def image_page(self, image_path: str | Path, title: str) -> None:
        self.title(title, size=13)
        image = ImageReader(str(Path(image_path).resolve()))
        image_width, image_height = image.getSize()
        available_width = self.width - 56
        available_height = self.height - 82
        scale = min(available_width / image_width, available_height / image_height)
        draw_width = image_width * scale
        draw_height = image_height * scale
        self.canvas.drawImage(
            image,
            (self.width - draw_width) / 2,
            31 + (available_height - draw_height) / 2,
            width=draw_width,
            height=draw_height,
            preserveAspectRatio=True,
            mask="auto",
        )
        self.finish_page()

    def paragraph(self, html: str, x: float, y: float, width: float, height: float,
                  *, font_size: float = 9.5, leading: float = 13.0,
                  light: bool = False) -> None:
        style = ParagraphStyle(
            "formal-report-paragraph",
            parent=self.styles["BodyText"],
            fontName="Helvetica",
            fontSize=font_size,
            leading=leading,
            textColor=colors.HexColor(self.theme.light_text if light else self.theme.navy),
        )
        paragraph = Paragraph(html, style)
        paragraph.wrapOn(self.canvas, width, height)
        paragraph.drawOn(self.canvas, x, y)

    def paragraph_top(self, html: str, x: float, top: float, width: float, height: float,
                      *, font_size: float = 9.5, leading: float = 13.0,
                      light: bool = False) -> None:
        style = ParagraphStyle(
            "formal-report-top-paragraph",
            parent=self.styles["BodyText"],
            fontName="Helvetica",
            fontSize=font_size,
            leading=leading,
            textColor=colors.HexColor(self.theme.light_text if light else self.theme.navy),
        )
        paragraph = Paragraph(html, style)
        _, used_height = paragraph.wrapOn(self.canvas, width, height)
        paragraph.drawOn(self.canvas, x, top - used_height)

    def table(self, rows: Sequence[Sequence[object]], *, x: float, y: float,
              column_widths: Sequence[float], row_heights: float | Sequence[float] = 23,
              font_size: float = 7.5) -> None:
        table = Table(rows, colWidths=list(column_widths), rowHeights=row_heights)
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(self.theme.blue)),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor(self.theme.pale)),
            ("TEXTCOLOR", (0, 1), (-1, -1), colors.HexColor(self.theme.navy)),
            ("GRID", (0, 0), (-1, -1), 0.45, colors.HexColor(self.theme.border)),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
            ("FONTSIZE", (0, 0), (-1, -1), font_size),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("ALIGN", (1, 0), (-1, -1), "CENTER"),
        ]))
        table.wrapOn(self.canvas, sum(column_widths), self.height)
        table.drawOn(self.canvas, x, y)

    def save(self) -> Path:
        self.canvas.save()
        return self.output


__all__ = ["FormalLandscapePDF", "FormalReportTheme"]
