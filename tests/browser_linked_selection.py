"""Manual browser regression: uv run --with playwright python tests/browser_linked_selection.py URL."""

import argparse
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from playwright.sync_api import TimeoutError as PlaywrightTimeoutError, expect, sync_playwright


OUTLINED_IDS = """el => {
    let ids = [];
    for (const trace of el.data) {
        if (!trace.meta?.observations) continue;
        const widths = trace.marker?.line?.width;
        if (widths === undefined) continue;
        for (let i = 0; i < trace.customdata.length; i++) {
            const width = typeof widths === 'number' ? widths : widths[i];
            if (width === 2.5) ids.push(trace.customdata[i][0]);
        }
    }
    return ids;
}"""


def row_checkbox(table, point_id):
    return table.locator("tr").filter(has_text=point_id).get_by_test_id("select-row-checkbox")


def wait_for_ids(page, plot, ids):
    try:
        page.wait_for_function(
            f"""([el, expected]) => {{
                const actual = [...new Set(({OUTLINED_IDS})(el))].sort();
                return JSON.stringify(actual) === JSON.stringify(expected.sort());
            }}""",
            arg=[plot.element_handle(), ids],
            timeout=30000,
        )
    except PlaywrightTimeoutError:
        print("Expected outlines:", ids, "actual:", plot.evaluate(OUTLINED_IDS))
        print("Checked rows:", page.locator("marimo-table").first.locator(
            'tr[data-state="selected"] [data-cell-id$="_point_id"]'
        ).all_text_contents())
        raise


def download(page, link):
    with page.expect_download() as event:
        link.click()
    return Path(event.value.path()).read_bytes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url", help="Running focused notebook URL, including its access token")
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    expect.set_options(timeout=30000)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(args=["--no-sandbox"])
        page = browser.new_page(viewport={"width": 1450, "height": 1100}, accept_downloads=True)
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(args.url)
        page.get_by_role("button", name="Load snapshot").click()
        plot = page.locator(".js-plotly-plot")
        page.wait_for_selector(".js-plotly-plot .scatterlayer", timeout=90000)
        original_count = int(re.search(
            r"(\d+) (?:regions|observations|points) shown", page.locator("body").inner_text()
        )[1])
        fixed_download = page.locator("a[download$='reconstruction.zip']")
        if not fixed_download.count():
            fixed_download = page.locator("a[download$='reconstruction.svg']")
        reconstruction = download(page, fixed_download)
        plot.evaluate("""el => {
            el.inspectionTest = true;
            const range = el._fullLayout.xaxis.range;
            const span = range[1]-range[0];
            window.Plotly.relayout(el, {'xaxis.range': [range[0]+span*.1,range[1]-span*.1]});
            window.Plotly.restyle(el, {visible:'legendonly'}, [1]);
        }""")
        expected_range = plot.evaluate("el=>el.layout.xaxis.range")
        table = page.locator("marimo-table").first
        table.scroll_into_view_if_needed()
        wanted = table.locator('[data-cell-id$="_point_id"]').first.inner_text().strip()
        first = row_checkbox(table, wanted)
        first.check()
        wait_for_ids(page, plot, [wanted])
        assert plot.evaluate("el=>el.inspectionTest"), "Chart was remounted"
        assert plot.evaluate("el=>el.layout.xaxis.range") == expected_range
        assert plot.evaluate("el=>el.data[1].visible") == "legendonly"
        # A shared region in the two RMSD panels must be outlined in both traces.
        assert plot.evaluate("""(el, wanted) => el.data.filter(t=>t.meta?.observations)
            .every(t=>t.customdata.every((r,i)=>r[0]!==wanted
                || t.marker.line.width[i]===2.5))""", wanted)
        assert download(page, fixed_download) == reconstruction
        plot.scroll_into_view_if_needed()
        if args.screenshot:
            args.screenshot.parent.mkdir(parents=True, exist_ok=True)
            plot.screenshot(path=str(args.screenshot))
        assert plot.evaluate("el=>el._context.toImageButtonOptions.format") == "svg"
        svg = download(page, plot.locator(".modebar-btn[data-title^='Download plot']"))
        ET.fromstring(svg)
        first.uncheck()
        wait_for_ids(page, plot, [])

        # Native box selection still filters the table and permits inspection inside it.
        page.reload()
        page.get_by_role("button", name="Load snapshot").click()
        page.wait_for_selector(".js-plotly-plot .scatterlayer", timeout=90000)
        plot.scroll_into_view_if_needed()
        plot.evaluate("""el => el.on('plotly_selected', event => {
            el.inspectionSelectedIds = [...new Set(event.points.map(p =>
                p.customdata?.[0] || el.data[p.curveNumber].customdata[p.pointIndex]?.[0]
            ).filter(Boolean))];
        })""")
        plot.locator(".modebar-btn[data-title='Box Select']").click()
        bounds = plot.evaluate("""el => {
            const l=el._fullLayout, b=el.getBoundingClientRect();
            return {x:b.x+l.xaxis._offset,y:b.y+l.yaxis._offset,
                    w:l.xaxis._length,h:l.yaxis._length};
        }""")
        page.mouse.move(bounds["x"] + bounds["w"]*.05, bounds["y"] + bounds["h"]*.2)
        page.mouse.down()
        page.mouse.move(bounds["x"] + bounds["w"]*.8, bounds["y"] + bounds["h"]*.98, steps=12)
        page.mouse.up()
        expect(page.get_by_test_id("select-row-checkbox").first).not_to_be_checked()
        page.wait_for_function(
            """el => el.inspectionSelectedIds?.length""", arg=plot.element_handle()
        )
        expected_ids = plot.evaluate("el => el.inspectionSelectedIds")
        expect(page.get_by_text(re.compile(
            rf"^{len(expected_ids)} (?:regions|observations|points) shown", re.MULTILINE
        ))).to_be_visible(timeout=30000)
        table.scroll_into_view_if_needed()
        expect(page.locator("marimo-table").first.locator(
            '[data-cell-id$="_point_id"]'
        ).first).to_have_text(re.compile(
            "^(?:" + "|".join(re.escape(point_id) for point_id in expected_ids) + ")$"
        ), timeout=30000)
        filtered_id = page.locator("marimo-table").first.locator(
            '[data-cell-id$="_point_id"]'
        ).first.inner_text().strip()
        first = row_checkbox(page.locator("marimo-table").first, filtered_id)
        first.check()
        assert filtered_id in expected_ids
        wait_for_ids(page, plot, [filtered_id])
        page.get_by_role("button", name="Clear plot selection", exact=True).click()
        expect(page.get_by_text(re.compile(
            rf"^{original_count} (?:regions|observations|points) shown", re.MULTILINE
        ))).to_be_visible(timeout=30000)
        wait_for_ids(page, plot, [])
        assert not errors, errors
        print(f"PASS: {original_count} observations; table→plot, clear, "
              f"preserved view, SVG, unchanged reconstruction, plot→table ({len(expected_ids)} selected)")
        browser.close()


if __name__ == "__main__":
    main()
