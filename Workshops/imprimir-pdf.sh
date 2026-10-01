#!/usr/bin/env bash
#
# Imprime el documento a PDF. No hay una fuente aparte para el PDF: es la misma
# página, así que las dos no pueden desincronizarse. Las reglas de impresión están
# al final del <style> del propio HTML.
#
#     ./imprimir-pdf.sh
#
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHROME="$(command -v google-chrome || command -v chromium || command -v chromium-browser || true)"

if [[ -z "$CHROME" ]]; then
    echo "No encuentro Chrome ni Chromium, que es lo que imprime." >&2
    exit 1
fi

"$CHROME" \
    --headless \
    --disable-gpu \
    --no-sandbox \
    --no-pdf-header-footer \
    --virtual-time-budget=15000 \
    --print-to-pdf="${DIR}/cuenta-earth-engine.pdf" \
    "file://${DIR}/cuenta-earth-engine.html" 2>/dev/null

# Chrome imprime lo que le den, incluida una página de error, y no dice nada.
if command -v pdfinfo >/dev/null; then
    pages="$(pdfinfo "${DIR}/cuenta-earth-engine.pdf" 2>/dev/null | awk '/^Pages:/ {print $2}')"
    if [[ "${pages:-0}" -lt 3 ]]; then
        echo "Ha salido con ${pages:-0} página(s): eso no es el documento." >&2
        exit 1
    fi
    printf 'cuenta-earth-engine.pdf  %s páginas  %s\n' "$pages" "$(du -h "${DIR}/cuenta-earth-engine.pdf" | cut -f1)"
else
    echo "cuenta-earth-engine.pdf listo."
fi
