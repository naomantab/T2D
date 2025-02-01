document.getElementById("tableFilter").addEventListener("keyup", function() {
    let filter = this.value.toLowerCase();
    let rows = document.querySelectorAll("#snp-table tbody tr");

    rows.forEach(row => {
      let rowText = row.textContent.toLowerCase(); // Get all the text of the row
      if (rowText.includes(filter)) { 
        row.style.display = ""; // Show the row if it matches the filter
      } else {
        row.style.display = "none"; // Hide the row if it doesn't match
      }
    });
  });