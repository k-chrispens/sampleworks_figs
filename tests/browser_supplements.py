"""Browser regression for 06–08: pass the directory-server URL and notebook filename."""

import argparse
from io import BytesIO
import json
import re
from zipfile import ZipFile

from playwright.sync_api import expect, sync_playwright

from browser_linked_selection import download, row_checkbox, wait_for_ids


SHADOW_PLOTS = """
function plots(root=document) {
    let found=[...root.querySelectorAll('.js-plotly-plot')];
    for (const el of root.querySelectorAll('*')) {
        if (el.shadowRoot) found.push(...plots(el.shadowRoot));
    }
    return found;
}
"""


def expand(page, key, opened=True):
    button = page.get_by_role("button", name=key, exact=True)
    if (button.get_attribute("aria-expanded") == "true") != opened:
        button.click()
    expect(button).to_have_attribute("aria-expanded", str(opened).lower())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("url")
    parser.add_argument("notebook", choices=[
        "06_fit_supplements.py", "07_training_tradeoffs.py", "08_sampling_budgets.py"])
    args = parser.parse_args()
    expected = {"06": 6, "07": 10, "08": 2}[args.notebook[:2]]
    first_key = {"06": "rscc", "07": "spirals", "08": "recycles-steps"}[args.notebook[:2]]
    expect.set_options(timeout=60000)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(args=["--no-sandbox"])
        page = browser.new_page(viewport={"width": 1500, "height": 1100}, accept_downloads=True)
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(args.url)
        page.wait_for_selector("a")
        href = page.locator("a").evaluate_all(
            "(es, name)=>es.find(e=>e.href.includes(name)).href", args.notebook)
        page.goto(href)
        expect(page.locator(".js-plotly-plot")).to_have_count(expected, timeout=150000)
        for plot in page.locator(".js-plotly-plot").all():
            page.wait_for_function("el=>el.data && el.layout.meta", arg=plot.element_handle())
            assert plot.evaluate("""el => {
                const ids=new Set(el.layout.meta.observations.map(p=>p.point_id));
                return el.data.every(t => !t.meta?.observations ||
                    t.customdata.every(p=>ids.has(p[0])));
            }""")
        plot = page.locator(".js-plotly-plot").first
        original_points = plot.evaluate("el=>el.layout.meta.observations")
        assert original_points
        expand(page, first_key)
        table = page.locator("marimo-table").filter(
            has=page.locator('[data-cell-id$="_point_id"]')).first
        table.scroll_into_view_if_needed()
        wanted = table.locator('[data-cell-id$="_point_id"]').first.inner_text().strip()
        plot.evaluate("""el=>{
            el.inspectionTest=true;
            const range=el._fullLayout.xaxis.range;
            window.Plotly.relayout(el, {'xaxis.range':[range[0]+.001,range[1]-.001]});
            window.Plotly.restyle(el, {visible:'legendonly'}, [el.data.length-1]);
        }""")
        expected_range = plot.evaluate("el=>el.layout.xaxis.range")
        row_checkbox(table, wanted).check()
        wait_for_ids(page, plot, [wanted])
        assert plot.evaluate("el=>el.inspectionTest")
        assert plot.evaluate("el=>el.layout.xaxis.range") == expected_range
        assert plot.evaluate("el=>el.data[el.data.length-1].visible") == "legendonly"
        # Restore the last trace so a selected point is visible in the camera export.
        plot.evaluate("el=>window.Plotly.restyle(el,{visible:true},[el.data.length-1])")
        plot.scroll_into_view_if_needed()
        svg = download(page, plot.locator(".modebar-btn[data-title^='Download plot']"))
        assert b"stroke-width: 2.5px" in svg, "Camera export lost inspection outlines"
        assert plot.evaluate("el=>el.layout.meta.observations") == original_points
        page.get_by_role("button", name="Build SVG + data bundle", exact=True).click()
        fixed = page.locator("a[download$='reconstruction.zip']")
        expect(fixed).to_be_visible(timeout=120000)
        with ZipFile(BytesIO(download(page, fixed))) as archive:
            recipe = json.loads(archive.read("recipe.json"))
            assert recipe["figure"]["layout"]["meta"]["observations"] == original_points
            assert b"stroke-width: 2.5px" not in archive.read("figure.svg")
            if args.notebook.startswith("06"):
                assert "highlighted-regions.csv" in archive.namelist()
            if args.notebook.startswith("08"):
                assert "contributors.csv" in archive.namelist()
        page.get_by_role("button", name="Clear plot selection", exact=True).first.click()
        wait_for_ids(page, plot, [])

        # Real browser box selection: summary/guide vertices must not become table records.
        plot.scroll_into_view_if_needed()
        plot.evaluate("""el=>el.on('plotly_selected',event=>{
            el.selectedTest=[...new Set(event.points.map(p=>p.customdata?.[0]).filter(Boolean))];
        })""")
        plot.locator(".modebar-btn[data-title='Box Select']").click()
        bounds = plot.evaluate("""el=>{
            const l=el._fullLayout,b=el.getBoundingClientRect();
            return {x:b.x+l.xaxis._offset,y:b.y+l.yaxis._offset,
                    w:l.xaxis._length,h:l.yaxis._length};
        }""")
        page.mouse.move(bounds["x"] + 2, bounds["y"] + 2)
        page.mouse.down()
        page.mouse.move(bounds["x"] + bounds["w"] - 2, bounds["y"] + bounds["h"] - 2, steps=12)
        page.mouse.up()
        page.wait_for_function("el=>el.selectedTest?.length", arg=plot.element_handle())
        selected = plot.evaluate("el=>el.selectedTest")
        print("Native box selected", len(selected), "IDs", flush=True)
        table.scroll_into_view_if_needed()
        expect(table.locator('[data-cell-id$="_point_id"]').first).to_have_text(
            re.compile("^(?:" + "|".join(re.escape(i) for i in selected) + ")$"))
        print("Filtered table ready; mounted plots:", page.locator(".js-plotly-plot").count(), flush=True)
        page.get_by_role("button", name="Clear plot selection", exact=True).first.scroll_into_view_if_needed()
        wait_for_ids(page, plot, [])
        expand(page, first_key, False)
        # Exercise every observation-bearing panel, not merely the first plot.
        for other in page.locator(".js-plotly-plot").all():
            key = other.evaluate("el=>el.layout.meta.inspection_key")
            print("Inspect", key, flush=True)
            if not other.evaluate("el=>el.layout.meta.observations.length"):
                continue
            expand(page, key)
            rows = page.locator("marimo-table").filter(
                has=page.locator('[data-cell-id$="_point_id"]:visible')).first
            rows.scroll_into_view_if_needed()
            identity = rows.locator('[data-cell-id$="_point_id"]').first.inner_text().strip()
            row_checkbox(rows, identity).check()
            wait_for_ids(page, other, [identity])
            row_checkbox(rows, identity).uncheck()
            wait_for_ids(page, other, [])
            expand(page, key, False)
        for candidate in page.locator(".js-plotly-plot").all():
            if not candidate.locator(".legendtoggle").count():
                continue
            group = candidate.evaluate("el=>el.data.find(t=>t.showlegend && t.legendgroup)?.legendgroup")
            if not group:
                continue
            candidate.scroll_into_view_if_needed()
            candidate.locator(".legendtoggle").first.click()
            page.wait_for_function("""([el,group])=>el.data.filter(t=>t.legendgroup===group)
                .every(t=>t.visible==='legendonly')""", arg=[candidate.element_handle(), group])
            candidate.locator(".legendtoggle").first.click()
            break

        if not args.notebook.startswith("08"):
            groups = page.locator("select").filter(
                has=page.locator('option[value="Protenix variants"]'))
            for label, allowed in [
                ("Protenix variants", ["protenix", "protenix-mini", "protenix-tiny"]),
                ("RF3 variants", ["rf3", "rf3-chiral-off"]),
            ]:
                groups.select_option(label=label)
                key = "classes" if args.notebook.startswith("06") else "calibration"
                page.wait_for_function("""([key,allowed])=>{""" + SHADOW_PLOTS + """
                    const el=plots()
                        .find(el=>el.layout?.meta?.inspection_key===key);
                    const rows=el?.layout.meta.observations;
                    return rows?.length && rows.every(row=>allowed.includes(row.model_name));
                }""", arg=[key, allowed], timeout=120000, polling=250)
            model = page.locator("select").filter(
                has=page.locator('option[value="Protenix-Mini"]'))
            model.select_option(label="Protenix-Mini")
            page.wait_for_function("""()=>{""" + SHADOW_PLOTS + """
                const rows=plots()[0]?.layout?.meta?.observations;
                return rows?.length && rows.every(p=>p.model_name==='protenix-mini');
            }""", timeout=120000, polling=250)
            expand(page, first_key)
            fresh_table = page.locator("marimo-table").filter(
                has=page.locator('[data-cell-id$="_point_id"]:visible')).first
            expect(fresh_table.locator('[data-cell-id$="_point_id"]').first).to_contain_text("protenix-mini")
            identity = fresh_table.locator('[data-cell-id$="_point_id"]').first.inner_text().strip()
            row_checkbox(fresh_table, identity).check()
            wait_for_ids(page, plot, [identity])
        assert not errors, errors
        print(f"PASS {args.notebook}: {expected} Plotly panels; {len(original_points)} first-panel rows; "
              "table→plot, preserved view, camera SVG, fixed bundle, native box selection")
        browser.close()


if __name__ == "__main__":
    main()
