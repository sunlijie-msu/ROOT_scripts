import requests
import re
import pandas as pd
from datetime import datetime

def list_recent_nds_mass_chains(years_back=5):
    """
    Fetches Nuclear Data Sheets (ISSN 0090-3752) using bulletproof cursor pagination.
    Extracts all 'A=value' mass chains safely avoiding any URL malformation errors.
    """
    current_year = datetime.now().year
    start_year = current_year - years_back + 1 
    
    base_url = "https://api.crossref.org/journals/0090-3752/works"
    headers = {'User-Agent': 'mailto:researcher@example.com'}
    
    items =[]
    cursor = "*"
    rows = 1000
    
    print(f"Scanning CrossRef database for complete NDS articles from {start_year} to {current_year}...")
    
    # 1. Fetch metadata using safe dictionary params to completely avoid 400 Bad Request errors
    while True:
        params = {
            'rows': rows,
            'cursor': cursor
        }
        
        try:
            response = requests.get(base_url, headers=headers, params=params, timeout=20)
            response.raise_for_status()
            data = response.json()
        except requests.exceptions.RequestException as e:
            print(f"Network error: {e}")
            break
            
        current_items = data.get('message', {}).get('items',[])
        if not current_items:
            break
            
        items.extend(current_items)
        print(f"  -> Fetched {len(items)} total journal articles so far...")
        
        next_cursor = data.get('message', {}).get('next-cursor')
        if not next_cursor or len(current_items) < rows:
            break
        cursor = next_cursor

    # 2. Advanced Regex matching single and combined chains (e.g. "A=240", "A=267, 271 and 275")
    pattern = re.compile(r'\bA\s*=\s*(\d+(?:(?:\s*,\s*|\s+and\s+|\s*&\s*)\d+)*)', re.IGNORECASE)
    results =[]
    
    for item in items:
        # Safely handle missing title arrays
        title_list = item.get('title') or[]
        if not title_list:
            continue
        title = title_list[0]
        
        match = pattern.search(title)
        if match:
            # Safely determine Year via strict multi-fallback logic
            year = None
            for date_field in ['published-print', 'published-online', 'issued', 'published']:
                date_info = item.get(date_field) or {}
                date_parts = date_info.get('date-parts') or []
                if date_parts and date_parts[0] and date_parts[0][0]:
                    year = date_parts[0][0]
                    break
            
            # Local year filtering 
            if year is None or not (start_year <= year <= current_year):
                continue

            # Normalize mass string (e.g., "267 and 271" -> "267,271")
            raw_mass = match.group(1)
            mass_a = re.sub(r'\s+(and|&)\s+', ',', raw_mass, flags=re.IGNORECASE)
            mass_a = re.sub(r'\s+', '', mass_a)
            
            # Extract Authors gracefully handling NoneTypes
            authors_list = []
            for author in (item.get('author') or[]):
                given = author.get('given', '').strip()
                family = author.get('family', '').strip()
                name = author.get('name', '').strip()
                
                if given and family:
                    authors_list.append(f"{given} {family}")
                elif family:
                    authors_list.append(family)
                elif name:
                    authors_list.append(name)
            author_str = "; ".join(authors_list) if authors_list else "No Authors Listed"
            
            # Precise Pages Extraction safely converting to string
            pages = str(item.get('page') or '')
            total_pages = item.get('page-count')
            
            if total_pages is None and pages:
                page_parts = re.split(r'[-–—]', pages)
                if len(page_parts) == 2 and page_parts[0].isdigit() and page_parts[1].isdigit():
                    total_pages = int(page_parts[1]) - int(page_parts[0]) + 1
                elif len(page_parts) == 1 and page_parts[0].isdigit():
                    total_pages = "N/A" 
            elif total_pages is None:
                total_pages = "N/A"
                
            results.append({
                "Mass A": mass_a,
                "Year": year,
                "Authors": author_str,
                "Pages String": pages,
                "Total Pages": total_pages,
                "Title": title
            })
                
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(by=["Year", "Mass A"], ascending=[False, False]).reset_index(drop=True)
    return df


if __name__ == "__main__":
    # ========================================================
    # CONFIGURATION: Define how many years back to extract here
    # ========================================================
    YEARS_TO_EXTRACT = 5
    
    df_mass_chains = list_recent_nds_mass_chains(years_back=YEARS_TO_EXTRACT)
    
    if df_mass_chains is not None and not df_mass_chains.empty:
        print("\n" + "="*100)
        print(df_mass_chains.to_string())
        print("="*100)
        print(f"\nTotal A=value articles successfully extracted: {len(df_mass_chains)}")
        
        csv_filename = f"NDS_Mass_Chains_{YEARS_TO_EXTRACT}Years.csv"
        df_mass_chains.to_csv(csv_filename, index=False)
        print(f"Data safely exported to: {csv_filename}\n")
    else:
        print("No articles found matching the criteria.")

