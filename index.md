

Code: `[slide link](/presentation.html)`

Result:

[slide link](/presentation.html)


### Embed

Or you could embed the xaringan slide into a webpage using iframe:

Code:

```
## CSS styles
<style>
.resp-container {
    position: relative;
    overflow: hidden;
    padding-top: 56.25%;
}

.testiframe {
    position: absolute;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    border: 0;
}
</style>


## html iframe to embed the slide.

<div class="resp-container">
    <iframe class="testiframe" src="https://diegotrindade.github.io/presentation/presentation.html">
      Fallback text here for unsupporting browsers, of which there are scant few.
    </iframe>
</div>

```

Result:

<style>
.resp-container {
    position: relative;
    overflow: hidden;
    padding-top: 56.25%;
}

.testiframe {
    position: absolute;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    border: 0;
}
</style>

<div class="resp-container">
    <iframe class="testiframe" src="https://diegotrindade.github.io/presentation/presentation.html">
      Fallback text here for unsupporting browsers, of which there are scant few.
    </iframe>
</div>