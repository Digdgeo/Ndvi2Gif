# Earth Engine Authentication

Before using Ndvi2Gif, you need to authenticate with Google Earth Engine.

## Creating an Earth Engine Account

1. Go to [Google Earth Engine signup](https://signup.earthengine.google.com/)
2. Sign in with your Google account
3. Fill out the registration form
4. Wait for approval (usually instant for `.edu` emails, or within 1-2 days)

```{note}
Earth Engine is **free for research, education, and non-profit use**!
```

## Authentication Methods

### Method 1: Interactive Authentication (Recommended for First Time)

```python
import ee

# Trigger authentication flow
ee.Authenticate()

# Initialize Earth Engine
ee.Initialize()
```

This will:
1. Open a browser window
2. Ask you to sign in with Google
3. Generate a verification code
4. Store credentials locally for future use

### Method 2: Notebook Authentication

In Jupyter notebooks or Colab:

```python
import ee

# Authenticate using notebook mode
ee.Authenticate(force=True)  # Use force=True to re-authenticate

# Initialize
ee.Initialize()
```

### Method 3: Using Service Accounts (Advanced)

For production environments or automated workflows:

```python
import ee

# Path to your service account key file
service_account = 'your-service-account@project.iam.gserviceaccount.com'
credentials = ee.ServiceAccountCredentials(service_account, 'path/to/key.json')

# Initialize with service account
ee.Initialize(credentials)
```

```{note}
Service accounts require setting up credentials in Google Cloud Platform. See [EE documentation](https://developers.google.com/earth-engine/guides/service_account) for details.
```

## Testing Your Authentication

Verify that authentication works:

```python
import ee

# Initialize
ee.Initialize()

# Test with a simple computation
point = ee.Geometry.Point([-3.7, 40.4])  # Madrid, Spain

# Grab the first Sentinel-2 scene covering the point in June 2023
image = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
         .filterBounds(point)
         .filterDate('2023-06-01', '2023-07-01')
         .first())

# Get a value
value = image.select('B8').reduceRegion(
    reducer=ee.Reducer.mean(),
    geometry=point,
    scale=10
).getInfo()

print(f"✓ Authentication successful! NIR value: {value['B8']}")
```

## Authentication Storage

Credentials are stored in:
- **Linux/Mac**: `~/.config/earthengine/`
- **Windows**: `C:\Users\YourName\.config\earthengine\`

You only need to authenticate once per machine.

## Common Authentication Issues

### Issue: "Permission denied" errors

Try re-authenticating:

```python
import ee
ee.Authenticate(force=True)
ee.Initialize()
```

### Issue: "Project not found"

Specify your project explicitly:

```python
import ee
ee.Initialize(project='your-project-id')
```

Find your project ID at [Google Cloud Console](https://console.cloud.google.com/).

### Issue: Browser doesn't open during authentication

Use the command-line authentication:

```bash
earthengine authenticate
```

Or manually copy the authentication URL from the error message.

### Issue: Token expired

Re-run the authentication:

```python
import ee
ee.Authenticate(force=True)
```

## Google Colab Specifics

In Google Colab, use:

```python
import ee
ee.Authenticate()
ee.Initialize(project='your-project-id')
```

Colab requires specifying a project ID.

## Best Practices

1. **Don't share credentials**: Never commit `.config/earthengine/` to version control
2. **Re-authenticate periodically**: Tokens can expire after long periods
3. **Use service accounts for production**: Better for automated scripts
4. **Test after authentication**: Always verify with a simple query

## Next Steps

Now that you're authenticated, let's create your first analysis in the [Quick Start](quick_start.md) guide.
