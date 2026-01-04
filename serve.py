from fastapi import FastAPI, Form
from fastapi.responses import HTMLResponse
import positionDigest

app = FastAPI()
@app.get("/", response_class=HTMLResponse)
async def read_form():
    return """

